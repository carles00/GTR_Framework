#include "renderer.h"

#include <algorithm> //sort

#include "camera.h"
#include "../gfx/gfx.h"
#include "../gfx/shader.h"
#include "../gfx/mesh.h"
#include "../gfx/texture.h"
#include "../gfx/fbo.h"
#include "../pipeline/prefab.h"
#include "../pipeline/material.h"
#include "../pipeline/animation.h"
#include "../utils/utils.h"
#include "../extra/hdre.h"
#include "../core/ui.h"

#include "scene.h"


using namespace SCN;

//some globals
GFX::Mesh sphere;

Renderer::Renderer(const char* shader_atlas_filename)
{
	render_wireframe = false;
	render_boundaries = false;
	show_shadowmaps = false;
	enable_render_priority = true;
	render_mode = eRenderMode::DEFERRED;
	enable_normal_map = true;
	enable_occ = true;
	enable_specular = false;
	enable_shadows = true;
	show_gbuffers = false;
	pbr_is_active = false;
	buffers_to_show[0] = 0;
	buffers_to_show[1] = 1;
	buffers_to_show[2] = 2;
	buffers_to_show[3] = 3;

	scene = nullptr;
	skybox_cubemap = nullptr;
	gbuffers_fbo = nullptr;
	illumination_fbo = nullptr;

	render_order = std::vector<RenderCall>();
	lights = std::vector<LightEntity*>();
	shadow_atlas_width = 4096;
	shadow_atlas_height = 2048;
	shadowmap_width = shadow_atlas_width / 4;
	shadowmap_height = shadow_atlas_height / 2;

	shadow_atlas_fbo = new GFX::FBO();
	shadow_atlas_fbo->setDepthOnly(shadow_atlas_width, shadow_atlas_height);
	shadow_atlas = shadow_atlas_fbo->depth_texture;

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
	sphere.uploadToVRAM();
}

//Sets up the light and render_order vectors, renders the shadowmaps
void Renderer::setupScene(Camera* camera)
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;
	
	//clear lights vector
	lights.clear();

	//parse all scene nodes and asign them to rendercall or litghts
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible)
			continue;

		if (ent->getType() == SCN::eEntityType::PREFAB) {
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
			if (pent->prefab)
				//if it is a prefab create the renderCalls
				walkEntities(&pent->root, camera);

		}else if (ent->getType() == SCN::eEntityType::LIGHT){
			lights.push_back((SCN::LightEntity*)ent);
		}
	}

	if(enable_shadows)
		generateShadowmaps();
}

//Walks over every prefab on the scene and creates the render calls
void Renderer::walkEntities(SCN::Node* node, Camera* camera) {
	if (!node->visible)
		return;

	//compute global matrix
	Matrix44 node_model = node->getGlobalMatrix(true);
	//does this node have a mesh? then we must render it
	if (node->mesh && node->material)
	{
		//compute the bounding box of the object in world space (by using the mesh bounding box transformed to world space)
		BoundingBox world_bounding = transformBoundingBox(node_model, node->mesh->box);

		//if bounding box is inside the camera frustum then the object is probably visible
		if (camera->testBoxInFrustum(world_bounding.center, world_bounding.halfsize))
		{
			//create the render calls
			createRenderCall(node_model, node->mesh, node->material, camera->eye);
		}
	}

	//iterate recursively with children
	for (int i = 0; i < node->children.size(); ++i)
		walkEntities(node->children[i], camera);
}

//Creates a render call for a given entity and stores it on the render_order vector
void SCN::Renderer::createRenderCall(Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, vec3 camera_pos) {
	RenderCall rc;
	Vector3f nodepos = model.getTranslation();
	rc.model = model;
	rc.mesh = mesh;
	rc.material = material;
	rc.distance_to_camera = nodepos.distance(camera_pos);
	render_order.push_back(rc);
}

//Generates shadowmaps for each light that casts a shadow, renders to shadow_atlas_fbo
void Renderer::generateShadowmaps() {
	GFX::startGPULabel("Shadowmaps");
	eRenderMode previous_mode = render_mode;
	render_mode = FLAT;
	Camera camera;

	int j = 0;
	int i = 0;
	for(auto light : lights)
	{
		//4 columns 2 rows
		if (i > 3) {
			j++;
			i = 0;
		}

		if (!light->cast_shadows)
			continue;

		if (light->light_type != eLightType::SPOT && light->light_type != eLightType::DIRECTIONAL)
			continue;


		vec3 pos = light->root.model.getTranslation();
		vec3 front = light->root.model.rotateVector(vec3(0, 0, -1));
		vec3 up = vec3(0, 1, 0);
		camera.lookAt(pos, pos + front, up);

		if (light->light_type == eLightType::SPOT) {
			camera.setPerspective(light->cone_info.y * 2, 1, light->near_distance, light->max_distance);
		}

		if (light->light_type == eLightType::DIRECTIONAL) {
			float halfarea = light->area / 2;
			camera.setOrthographic(-halfarea, halfarea, halfarea * 1, -halfarea * 1, 0.1, light->max_distance);
		}
		
		//TODO check if light inside camera

		vec4 region = vec4(
			(shadowmap_width*i)/ shadow_atlas_width,
			(shadowmap_height*j)/ shadow_atlas_height,
			shadowmap_width / shadow_atlas_width, 
			shadowmap_height / shadow_atlas_height
		);
		light->shadowmap_region = region;

		shadow_atlas_fbo->bind();

		glViewport(region.x * shadow_atlas_width, region.y * shadow_atlas_height, region.z * shadow_atlas_width, region.w * shadow_atlas_height);
		glScissor(region.x * shadow_atlas_width, region.y * shadow_atlas_height, region.z * shadow_atlas_width, region.w * shadow_atlas_height);
		glEnable(GL_SCISSOR_TEST);

		renderFrame(&camera);
		
		glDisable(GL_SCISSOR_TEST);
		shadow_atlas_fbo->unbind();

		light->shadow_viewproj = camera.viewprojection_matrix;
		i++;
	}
	render_mode = previous_mode;
	GFX::endGPULabel();
}

//Render scene pipeline
void Renderer::renderScene(SCN::Scene* scene, Camera* camera)
{
	this->scene = scene;
	
	setupScene(camera);

	//render entities
	renderFrame(camera);

	//clear render_order after all rendering is finished
	render_order.clear();

	if (show_shadowmaps)
		showShadowmaps();

}

void Renderer::renderFrame(Camera* camera)
{
	if (render_mode == eRenderMode::DEFERRED)
		renderDeferred(camera);
	else
		renderForward(camera);
}



//Calls renderMesh for each rendercall stored in render_order
void SCN::Renderer::renderForward(Camera* camera) {
	

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	//set the camera as default (used by some functions in the framework)
	camera->enable();

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if (skybox_cubemap && render_mode != eRenderMode::FLAT)
		renderSkybox(skybox_cubemap);

	//order render_order for priority rendering
	if(enable_render_priority)
		std::sort(render_order.begin(), render_order.end(), RenderCall::render_call_sort);
	
	//call renderMeshWithMaterial depending on rendermode
	renderSceneNodes(camera);
}

void Renderer::renderSceneNodes(Camera* camera) {

	for (int i = 0; i < render_order.size(); i++) {
		switch (render_mode)
		{
		case FLAT:
			renderMeshWithMaterialFlat(&render_order[i], camera);
			break;
		case SCN::TEXTURED:
			renderMeshWithMaterial(&render_order[i], camera);
			break;
		case SCN::LIGHTS_MULTIPASS:
		case SCN::LIGHTS_SINGLEPASS:
			renderMeshWithMaterialLight(&render_order[i], camera);
			break;
		case SCN::DEFERRED:
			renderMeshWithMaterialGBuffers(&render_order[i], camera);
			break;
		default:
			break;
		}
	}
}


void Renderer::renderSkybox(GFX::Texture* cubemap)
{
	Camera* camera = Camera::current;

	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	GFX::Shader* shader = GFX::Shader::Get("skybox");
	if (!shader)
		return;
	shader->enable();

	Matrix44 m;
	m.setTranslation(camera->eye.x, camera->eye.y, camera->eye.z);
	m.scale(10, 10, 10);
	shader->setUniform("u_model", m);
	cameraToShader(camera, shader);
	shader->setUniform("u_texture", cubemap, 0);
	sphere.render(GL_TRIANGLES);
	shader->disable();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_DEPTH_TEST);
}

//renders a mesh given its transform and material
void Renderer::renderMeshWithMaterial(RenderCall* rc, Camera* camera)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material )
		return;
    assert(glGetError() == GL_NO_ERROR);

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* white = GFX::Texture::getWhiteTexture(); //a 1x1 white texture
	
	GFX::Texture* albedo_texture = rc->material->textures[SCN::eTextureChannel::ALBEDO].texture;

	//select the blending
	if (rc->material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if(rc->material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
    assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("texture");

    assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", rc->model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t );

	shader->setUniform("u_color", rc->material->color);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white, 0);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", rc->material->alpha_mode == SCN::eAlphaMode::MASK ? rc->material->alpha_cutoff : 0.001f);

	if (render_wireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	//do the draw call that renders the mesh into the screen
	rc->mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

//renders Mesh with flat shader
void Renderer::renderMeshWithMaterialFlat(RenderCall* rc, Camera* camera)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;

	//select the blending
	if (rc->material->alpha_mode == SCN::eAlphaMode::BLEND)
		return;
	glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (rc->material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("flat");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", rc->model);
	cameraToShader(camera, shader);

	//do the draw call that renders the mesh into the screen
	rc->mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();
}

//Render mesh with lights: singlepass and multipas
void Renderer::renderMeshWithMaterialLight(RenderCall* rc, Camera* camera)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* white = GFX::Texture::getWhiteTexture(); //a 1x1 white texture

	//get textures from material
	GFX::Texture* albedo_texture = rc->material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = rc->material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	GFX::Texture* metalic_texture = rc->material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	GFX::Texture* normal_map = rc->material->textures[SCN::eTextureChannel::NORMALMAP].texture;

	//select the blending
	if (rc->material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (rc->material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	if (render_mode == eRenderMode::LIGHTS_MULTIPASS) {
		shader = GFX::Shader::Get("light_multipass");
	}
	else if (render_mode == eRenderMode::LIGHTS_SINGLEPASS) {
		shader = GFX::Shader::Get("light_singlepass");
	}
	

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", rc->model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);
	//pass all textures and color to shader
	shader->setUniform("u_color", rc->material->color);
	shader->setUniform("u_emissive_factor", rc->material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white, 1);
	shader->setUniform("u_ambient", scene->ambient_light);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", rc->material->alpha_mode == SCN::eAlphaMode::MASK ? rc->material->alpha_cutoff : 0.001f);
	
	//set texture flags
	vec4 texture_flags = vec4(0,0,0,0); //normal, occlusion,specular | if 0 there is no texture 

	//----normal map
	if (normal_map && enable_normal_map) {
		texture_flags.x = 1;
		shader->setUniform("u_normal_map",normal_map,2);
	}
	//----occ & specular
	if (metalic_texture && (enable_occ|| enable_specular || pbr_is_active)) {
		shader->setUniform("u_metalic_roughness", vec2(rc->material->metallic_factor, rc->material->roughness_factor));
		shader->setUniform("u_view_pos", camera->eye);
		if(enable_occ)
			texture_flags.y = 1;

		if (enable_specular) 
			//disabled
			texture_flags.z = 0;
		
		if (pbr_is_active) 
			texture_flags.w = 1;
		

		shader->setUniform("u_occ_met_rough_texture", metalic_texture,3);
	}

	shader->setUniform("u_texture_flags", texture_flags);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//allow rendering at the same depth
	glDepthFunc(GL_LEQUAL);

	//select single or multipass
	if (render_mode == eRenderMode::LIGHTS_MULTIPASS) {
		multiPass(rc, shader);
	}
	else if (render_mode == eRenderMode::LIGHTS_SINGLEPASS) {
		singlePass(rc, shader);
	}
	

	glDepthFunc(GL_LESS);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

//Render function for multipass shader
void Renderer::multiPass(RenderCall* rc, GFX::Shader* shader) {
	// TODO BUG when there is no Directional light and no spotlight, everything that is not touched by a spot light does not rendered, but it should 
	if (lights.size() == 0) {
		shader->setUniform("u_light_info", vec4(int(eLightType::NO_LIGHT), 0, 0, 0));
		rc->mesh->render(GL_TRIANGLES);
		return;
	}

	for (int i = 0; i < lights.size(); i++)
	{
		//if distance from light is too big, continue
		LightEntity* light = lights[i];
		eLightType type = light->light_type;
		//only check spot and point because directional and ambient are global
		if (type == eLightType::POINT || type == eLightType::SPOT) {
			BoundingBox world_bounding = transformBoundingBox(rc->model, rc->mesh->box);
			if (cullLights(light, world_bounding))
				continue;
		}

		lightToShader(light, shader);

		//do the draw call that renders the mesh into the screen
		rc->mesh->render(GL_TRIANGLES);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		shader->setUniform("u_ambient", vec3(0.0));
		shader->setUniform("u_emissive_factor", vec3(0.0));
	}
}

//Render function for singlepass shader
void Renderer::singlePass(RenderCall* rc, GFX::Shader* shader) {
	if (lights.size() == 0) {
		shader->setUniform("u_num_lights", 0);
		rc->mesh->render(GL_TRIANGLES);
		return;
	}
	//lights
	vec4 lights_info[MAX_LIGHTS];
	vec3 lights_fronts[MAX_LIGHTS];
	vec3 lights_positions[MAX_LIGHTS];
	vec3 lights_colors[MAX_LIGHTS];
	vec2 lights_cones[MAX_LIGHTS];
	//shadows
	vec2 shadow_params[MAX_LIGHTS];
	mat4 shadow_viewproj[MAX_LIGHTS];
	vec4 shadowmap_region[MAX_LIGHTS];

	for (int i = 0; i < MAX_LIGHTS; i++)
	{
		if (i < lights.size()) {
			//lights
			LightEntity* light = lights[i];
			lights_info[i] = vec4((int)light->light_type, light->near_distance, light->max_distance, 0);
			lights_fronts[i] = light->root.model.rotateVector(vec3(0, 0, 1));
			lights_positions[i] = light->root.model.getTranslation();
			lights_colors[i] = light->color * light->intensity;
			if (light->light_type == eLightType::SPOT)
				lights_cones[i] = vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD));
			else
				lights_cones[i] = vec2(0, 0);
			//shadows
			shadow_params[i] = vec2(light->cast_shadows && enable_shadows ? 1 : 0, light->shadow_bias);
			if (light->cast_shadows && enable_shadows) {
				shadow_viewproj[i] = light->shadow_viewproj;
				shadowmap_region[i] = light->shadowmap_region;
			}
		}
	}
	int size = lights.size();
	//lights
	shader->setUniform4Array("u_light_info", (float*)&lights_info, size);
	shader->setUniform3Array("u_light_front", (float*)&lights_fronts, size);
	shader->setUniform3Array("u_light_position", (float*)&lights_positions, size);
	shader->setUniform3Array("u_light_color", (float*)&lights_colors, size);
	shader->setUniform2Array("u_light_cone", (float*)&lights_cones, size);
	shader->setUniform1("u_num_lights", size);

	//shadows
	shader->setUniform2Array("u_shadow_params", (float*)&shadow_params, size);
	if (enable_shadows) {
		shader->setUniform4Array("u_shadow_region", (float*)&shadowmap_region, size);
		shader->setMatrix44Array("u_shadow_viewproj", shadow_viewproj, size);
		shader->setTexture("u_shadowmap", shadow_atlas, 8);
	}
	//do the draw call that renders the mesh into the screen

	rc->mesh->render(GL_TRIANGLES);
}

void SCN::Renderer::renderDeferred(Camera* camera) {

	vec2 size = CORE::getWindowSize();
	//generate the gbuffers
	if (!gbuffers_fbo) {

		gbuffers_fbo = new GFX::FBO();
		gbuffers_fbo->create(size.x, size.y, 4, GL_RGBA, GL_UNSIGNED_BYTE, true);

	}
	gbuffers_fbo->bind();
	//Renderingi inside the gBuffer
	{
		/*gbuffers_fbo->enableBuffers(true, false, false, false);*/
		gbuffers_fbo->enableAllBuffers();
		glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		renderSceneNodes(camera);
	}
	gbuffers_fbo->unbind();
	
	
	camera->enable();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(skybox_cubemap)
		renderSkybox(skybox_cubemap);
	//glDisable(GL_DEPTH_TEST);
	GFX::Mesh* quad = GFX::Mesh::getQuad();
	
	GFX::Shader* quad_shader = GFX::Shader::Get("deferred_global");
	quad_shader->enable();

	
	quad_shader->setTexture("u_albedo_texture", gbuffers_fbo->color_textures[0], 0);
	quad_shader->setTexture("u_normal_texture", gbuffers_fbo->color_textures[1], 1);
	quad_shader->setTexture("u_extra_texture", gbuffers_fbo->color_textures[2], 2);
	quad_shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);

	quad_shader->setUniform("u_ambient_light", scene->ambient_light);

	quad->render(GL_TRIANGLES);

	GFX::Shader* light_shader = GFX::Shader::Get("deferred_light");
	light_shader->enable();
	light_shader->setTexture("u_albedo_texture", gbuffers_fbo->color_textures[0], 0);
	light_shader->setTexture("u_normal_texture", gbuffers_fbo->color_textures[1], 1);
	light_shader->setTexture("u_extra_texture", gbuffers_fbo->color_textures[2], 2);
	light_shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);

	if (gbuffers_fbo->color_textures[3] && pbr_is_active)
	{
		light_shader->setTexture("u_metalic_roughness", gbuffers_fbo->color_textures[3], 4);
		light_shader->setUniform("u_pbr_state", 1.0f);
	}else
		light_shader->setUniform("u_pbr_state", 0.0f);

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	vec3 ambient = scene->ambient_light;
	for (auto light : lights) {
		
		light_shader->setUniform("u_ivp", camera->inverse_viewprojection_matrix);
		light_shader->setUniform("u_iRes", vec2(1.0 / size.x, 1.0 / size.y));
		light_shader->setUniform("u_ambient_light", ambient);
		light_shader->setUniform("u_eye", camera->eye);
		lightToShader(light, light_shader);
		quad->render(GL_TRIANGLES);
		ambient = vec3(0.0f); 
		
	}
	glDisable(GL_BLEND);

	//TODO: Create sphere for  every light
	// such that renders only the parts INSIDE the mesh
	/*glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	GFX::Shader* shader = GFX::Shader::Get("flat");
	shader->enable();
	Matrix44 model;
	model.scale(10, 10, 10);
	shader->setUniform("u_model", model);
	shader->setUniform("u_color", vec4(1.0,0.0,0.0,1.0));
	cameraToShader(camera, shader);
	sphere.render(GL_TRIANGLES);*/

	if (show_gbuffers)
	{
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
		float halfWidth, halfHeight;
		halfWidth = size.x / 2;
		halfHeight = size.y / 2;
		glViewport(0, halfHeight, halfWidth, halfHeight);
		if (buffers_to_show[0] < 3)
			gbuffers_fbo->color_textures[buffers_to_show[0]]->toViewport();
		else if (buffers_to_show[0] > 3)
		{
			GFX::Shader* texture_shader = nullptr;
			switch (buffers_to_show[0]) {
				case 4:
					texture_shader = GFX::Shader::getDefaultShader("screen_channel_r");
					texture_shader->enable();
					gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
					break;
				case 5:
					texture_shader = GFX::Shader::getDefaultShader("screen_channel_g");
					texture_shader->enable();
					gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
					break;
				case 6:
					texture_shader = GFX::Shader::getDefaultShader("screen_channel_b");
					texture_shader->enable();
					gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
					break;
			}
			
		}
		else
		{
			GFX::Shader* deph_shader = GFX::Shader::getDefaultShader("linear_depth");
			deph_shader->enable();
			gbuffers_fbo->depth_texture->toViewport(deph_shader);
			deph_shader->disable();
		}
		glViewport(halfWidth, halfHeight, halfWidth, halfHeight);
		if (buffers_to_show[1] < 3)
			gbuffers_fbo->color_textures[buffers_to_show[1]]->toViewport();
		else if (buffers_to_show[1] > 3)
		{
			GFX::Shader* texture_shader = nullptr;
			switch (buffers_to_show[1]) {
			case 4:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_r");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 5:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_g");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 6:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_b");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			}
			texture_shader->disable();
		}
		else
		{
			GFX::Shader* deph_shader = GFX::Shader::getDefaultShader("linear_depth");
			deph_shader->enable();
			gbuffers_fbo->depth_texture->toViewport(deph_shader);
			deph_shader->disable();
		}
		glViewport(0, 0, halfWidth, halfHeight);
		if (buffers_to_show[2] < 3)
			gbuffers_fbo->color_textures[buffers_to_show[2]]->toViewport();
		else if (buffers_to_show[2] > 3)
		{
			GFX::Shader* texture_shader = nullptr;
			switch (buffers_to_show[2]) {
			case 4:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_r");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 5:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_g");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 6:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_b");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			}
			texture_shader->disable();
		}
		else
		{
			GFX::Shader* deph_shader = GFX::Shader::getDefaultShader("linear_depth");
			deph_shader->enable();
			gbuffers_fbo->depth_texture->toViewport(deph_shader);
			deph_shader->disable();
		}
		glViewport(halfWidth, 0, halfWidth, halfHeight);
		/*gbuffers_fbo->color_textures[3]->toViewport();
		glViewport(halfWidth / 2, halfHeight / 2, halfWidth, halfHeight);*/
		GFX::Shader* deph_shader = GFX::Shader::getDefaultShader("linear_depth");
		deph_shader->enable();
		deph_shader->setUniform("u_camera_nearfar", vec2(camera->near_plane, camera->far_plane));
		if (buffers_to_show[3] < 3)
			gbuffers_fbo->color_textures[buffers_to_show[2]]->toViewport();
		else if (buffers_to_show[3] > 3)
		{
			GFX::Shader* texture_shader = nullptr;
			switch (buffers_to_show[3]) {
			case 4:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_r");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 5:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_g");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			case 6:
				texture_shader = GFX::Shader::getDefaultShader("screen_channel_b");
				texture_shader->enable();
				gbuffers_fbo->color_textures[3]->toViewport(texture_shader);
				break;
			}
			texture_shader->disable();
		}
		else
		{
			GFX::Shader* deph_shader = GFX::Shader::getDefaultShader("linear_depth");
			deph_shader->enable();
			gbuffers_fbo->depth_texture->toViewport(deph_shader);
			deph_shader->disable();
		}
		glViewport(0, 0, size.x, size.y);
		deph_shader->disable();
	}
	//compute the illumination


}
void Renderer::renderMeshWithMaterialGBuffers(RenderCall* rc, Camera* camera)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (rc->material->alpha_mode == eAlphaMode::BLEND)
		return;

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* white = GFX::Texture::getWhiteTexture(); //a 1x1 white texture

	GFX::Texture* albedo_texture = rc->material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = rc->material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	GFX::Texture* normal_texture = rc->material->textures[SCN::eTextureChannel::NORMALMAP].texture;
	GFX::Texture* metalness_texture = rc->material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;

	//select the blending
	glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (rc->material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("gbuffers");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", rc->model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);
	//TODO: normalmaps with flags
	vec4 texture_flags = vec4(0, 0, 0, 0); //normal, occlusion, specular
	
	shader->setUniform("u_color", rc->material->color);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white, 1);
	shader->setUniform("u_emissive_factor", rc->material->emissive_factor);
	int textureSlot = 2;
	if (normal_texture && enable_normal_map) {
		texture_flags.x = 1;
		shader->setUniform("u_normal_texture", normal_texture, textureSlot++);
	}
	if (metalness_texture)
		shader->setUniform("u_metalic_roughness", metalness_texture, textureSlot++);

	shader->setUniform("u_alpha_cutoff", (float) (rc->material->alpha_mode == SCN::eAlphaMode::MASK ? rc->material->alpha_cutoff : 0.001f));
	
	shader->setUniform("u_texture_flags", texture_flags);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//do the draw call that renders the mesh into the screen
	rc->mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void SCN::Renderer::cameraToShader(Camera* camera, GFX::Shader* shader)
{
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix );
	shader->setUniform("u_camera_position", camera->eye);
}
void SCN::Renderer::lightToShader(LightEntity* light, GFX::Shader* shader)
{
	shader->setUniform("u_light_info", vec4((int)light->light_type, light->near_distance, light->max_distance, 0));
	shader->setUniform("u_light_front", light->root.model.rotateVector(vec3(0, 0, 1)));
	shader->setUniform("u_light_position", light->root.model.getTranslation());
	shader->setUniform("u_light_color", light->color * light->intensity);
	if (light->light_type == eLightType::SPOT)
		shader->setUniform("u_light_cone", vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD)));


	shader->setUniform("u_shadow_params", vec2(light->cast_shadows && enable_shadows ? 1 : 0, light->shadow_bias));
	if (light->cast_shadows && enable_shadows) {
		shader->setTexture("u_shadowmap", shadow_atlas, 8);
		shader->setUniform("u_shadow_viewproj", light->shadow_viewproj);
		shader->setUniform("u_shadow_region", light->shadowmap_region);
	}
}

void Renderer::showShadowmaps() {
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	if (!enable_shadows)
		return;

	GFX::Shader* shader = GFX::Shader::getDefaultShader("linear_depth");
	shader->enable();
	shader->setUniform("u_camera_nearfar", vec2( 0.1, 1000));
	glViewport(300, 100, 512, 256);
	shadow_atlas->toViewport( shader);
	shader->disable();
		
	
	vec2 size = CORE::getWindowSize();
	glViewport(0, 0, size.x, size.y);
}
bool Renderer::cullLights(LightEntity* light, BoundingBox bb) {

	eLightType type = light->light_type;
	if (type == eLightType::POINT )
		return !BoundingBoxSphereOverlap(bb, light->root.model.getTranslation(), light->max_distance);
	if (type == eLightType::SPOT) {
		return !spotLightAABB(light, bb);
	}
}

bool Renderer::spotLightAABB(LightEntity* light, BoundingBox bb) {
	float sphereRadius = max(bb.halfsize.x * 2.0f, bb.halfsize.z * 2.0f);
	vec3 v = bb.center - light->root.model.getTranslation();
	float lenSq = dot(v, v);
	float v1Len = dot(v, light->root.model.rotateVector(vec3(0, 0, -1)));
	
	float distanceClosestPoint = cos(light->cone_info.y * DEG2RAD) * sqrt(lenSq - v1Len * v1Len) - v1Len * sin(light->cone_info.y * DEG2RAD);
	bool angleCull = distanceClosestPoint > sphereRadius;
	bool frontCull = v1Len > sphereRadius + light->max_distance;
	bool backCull = v1Len < -sphereRadius;

	return !(angleCull || frontCull || backCull);
}

//Switch widget extracted from: https://github.com/ocornut/imgui/issues/1537#issuecomment-780262461
void ToggleButton(const char* label, bool* v)
{
	ImGui::Text(label);
	ImGui::SameLine();
	ImVec4* colors = ImGui::GetStyle().Colors;
	ImVec2 p = ImGui::GetCursorScreenPos();
	ImDrawList* draw_list = ImGui::GetWindowDrawList();
	const char* str_id = "";

	float height = ImGui::GetFrameHeight();
	float width = height * 1.55f;
	float radius = height * 0.50f;

	ImGui::InvisibleButton(str_id, ImVec2(width, height));
	if (ImGui::IsItemClicked()) *v = !*v;
	ImGuiContext& gg = *GImGui;
	float ANIM_SPEED = 0.085f;
	if (gg.LastActiveId == gg.CurrentWindow->GetID(str_id))// && g.LastActiveIdTimer < ANIM_SPEED)
		float t_anim = ImSaturate(gg.LastActiveIdTimer / ANIM_SPEED);
	if (ImGui::IsItemHovered())
		draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), ImGui::GetColorU32(*v ? colors[ImGuiCol_ButtonActive] : ImVec4(0.78f, 0.78f, 0.78f, 1.0f)), height * 0.5f);
	else
		draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), ImGui::GetColorU32(*v ? colors[ImGuiCol_Button] : ImVec4(0.85f, 0.85f, 0.85f, 1.0f)), height * 0.50f);
	draw_list->AddCircleFilled(ImVec2(p.x + radius + (*v ? 1 : 0) * (width - radius * 2.0f), p.y + radius), radius - 1.5f, IM_COL32(255, 255, 255, 255));
}
#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ToggleButton("Wireframe", &render_wireframe);
	ToggleButton("Boundaries", &render_boundaries);

	//add here your stuff
	ToggleButton("Render Priority", &enable_render_priority);
	ToggleButton("Normal Map", &enable_normal_map);
	ToggleButton("Occlusion", &enable_occ);
	ToggleButton("Specular", &enable_specular);
	ToggleButton("Shadows", &enable_shadows);
	ToggleButton("Show Shadow Atlas", &show_shadowmaps);
	ToggleButton("Activate PBR", &pbr_is_active);
	ImGui::Combo("Render Mode",(int*) &render_mode, "FLAT\0TEXTURED\0LIGHTS_MULTIPASS\0LIGHTS_SINGLEPASS\0DEFERRED\0", 5);
	if (render_mode == eRenderMode::DEFERRED)
		ToggleButton("Show all buffers", &show_gbuffers);
	if (show_gbuffers) {
		ImGui::Text("Select buffers");
		
		int buffer1 = buffers_to_show[0];
		int buffer2 = buffers_to_show[1];
		int buffer3 = buffers_to_show[2];
		int buffer4 = buffers_to_show[3];


		ImGui::Combo("Upper Left", (int*) &buffer1, "Albedo\0Normalmap\0Emissive\0Depth\0Oclussion\0Metallic\0Roughness", 7);
		ImGui::Combo("Upper Right", (int*) &buffer2, "Albedo\0Normalmap\0Emissive\0Depth\0Oclussion\0Metallic\0Roughness", 7);
		ImGui::Combo("Lower Left", (int*) &buffer3, "Albedo\0Normalmap\0Emissive\0Depth\0Oclussion\0Metallic\0Roughness", 7);
		ImGui::Combo("Upper Right", (int*) &buffer4, "Albedo\0Normalmap\0Emissive\0Depth\0Oclussion\0Metallic\0Roughness", 7);

		buffers_to_show[0] = buffer1;
		buffers_to_show[1] = buffer2;
		buffers_to_show[2] = buffer3;
		buffers_to_show[3] = buffer4;

		
	}
}

#else
void Renderer::showUI() {}
#endif