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
	enable_render_priority = true;
	render_mode = eRenderMode::LIGHTS_MULTIPASS;
	enable_normal_map = true;
	enable_occ = true;
	enable_specular = false;

	scene = nullptr;
	skybox_cubemap = nullptr;
	render_order = std::vector<RenderCall>();
	lights = std::vector<LightEntity*>();

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
}

void Renderer::setupScene()
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
				walkEntities(&pent->root, Camera::current);

		}else if (ent->getType() == SCN::eEntityType::LIGHT){
			lights.push_back((SCN::LightEntity*)ent);
		}
	}
}

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
			createRenderCall(node_model, node->mesh, node->material);
		}
	}

	//iterate recursively with children
	for (int i = 0; i < node->children.size(); ++i)
		walkEntities(node->children[i], camera);
}


void SCN::Renderer::createRenderCall(Matrix44 model, GFX::Mesh* mesh, SCN::Material* material) {
	RenderCall rc;
	Vector3f nodepos = model.getTranslation();
	rc.model = model;
	rc.mesh = mesh;
	rc.material = material;
	rc.distance_to_camera = nodepos.distance(Camera::current->eye);
	render_order.push_back(rc);
}

void Renderer::renderScene(SCN::Scene* scene, Camera* camera)
{
	this->scene = scene;
	
	//set the camera as default (used by some functions in the framework)
	camera->enable();
	
	setupScene();

	//render entities
	renderFrame();

}

void SCN::Renderer::renderFrame() {
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if (skybox_cubemap)
		renderSkybox(skybox_cubemap);

	//order render_order
	if(enable_render_priority)
		std::sort(render_order.begin(), render_order.end(), RenderCall::render_call_sort);
	
	//call renderMeshWithMaterial
	for (int i = 0; i < render_order.size(); i++) {
		switch (render_mode)
		{
			case SCN::FLAT:
				renderMeshWithMaterial(&render_order[i]);
				break;
			case SCN::LIGHTS_MULTIPASS:
			case SCN::LIGHTS_SINGLEPASS:
				renderMeshWithMaterialLight(&render_order[i]);
				break;
			default:
				break;
		}
	}
	//once done empty render_order
	render_order.clear();
}

void Renderer::multiPass(RenderCall* rc, GFX::Shader* shader) {

	if (lights.size() == 0) {
		shader->setUniform("u_light_info", vec4(int(eLightType::NO_LIGHT), 0, 0, 0));
		rc->mesh->render(GL_TRIANGLES);
		return;
	}
	
	for (int i = 0; i < lights.size(); i++)
	{
		//if distance from light is too big, continue
		LightEntity* light = lights[i];
		//only check spot and point because directional and ambient are global
		if (light->light_type == eLightType::SPOT || light->light_type == eLightType::POINT) {
			BoundingBox world_bounding = transformBoundingBox(rc->model, rc->mesh->box);
			if (!BoundingBoxSphereOverlap(world_bounding, light->root.model.getTranslation(), light->max_distance)) {
				continue;
			}
		}

		shader->setUniform("u_light_info", vec4((int)light->light_type, light->near_distance, light->max_distance, 0));
		shader->setUniform("u_light_front", light->root.model.rotateVector(vec3(0, 0, 1)));
		shader->setUniform("u_light_position", light->root.model.getTranslation());
		shader->setUniform("u_light_color", light->color * light->intensity);
		if (light->light_type == eLightType::SPOT)
			shader->setUniform("u_light_cone", vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD)));

		//do the draw call that renders the mesh into the screen
		rc->mesh->render(GL_TRIANGLES);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		shader->setUniform("u_ambient", vec3(0.0));
		shader->setUniform("u_emissive_factor", vec3(0.0));
	}
}

void Renderer::singlePass(RenderCall* rc, GFX::Shader* shader) {
	if (lights.size() == 0) {
		shader->setUniform("u_num_lights", 0);
		rc->mesh->render(GL_TRIANGLES);
		return;
	}

	vec4 lights_info[MAX_LIGHTS];
	vec3 lights_fronts[MAX_LIGHTS];
	vec3 lights_positions[MAX_LIGHTS];
	vec3 lights_colors[MAX_LIGHTS];
	vec2 lights_cones[MAX_LIGHTS];

	for (int  i = 0; i < MAX_LIGHTS; i++)
	{
		if (i < lights.size()) {
			LightEntity* light = lights[i];
			lights_info[i] = vec4((int)light->light_type, light->near_distance, light->max_distance, 0);
			lights_fronts[i] = light->root.model.rotateVector(vec3(0, 0, 1));
			lights_positions[i] = light->root.model.getTranslation();
			lights_colors[i] = light->color * light->intensity;
			if (light->light_type == eLightType::SPOT)
				lights_cones[i] = vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD));
			else
				lights_cones[i] = vec2(0, 0);
		}
	}
	int size = lights.size();

	shader->setUniform4Array("u_light_info", (float *) &lights_info,size);
	shader->setUniform3Array("u_light_front", (float*) &lights_fronts, size);
	shader->setUniform3Array("u_light_position", (float*)&lights_positions, size);
	shader->setUniform3Array("u_light_color", (float*)&lights_colors, size);
	shader->setUniform2Array("u_light_cone", (float*)&lights_cones, size);

	shader->setUniform1("u_num_lights", size);

	//do the draw call that renders the mesh into the screen

	rc->mesh->render(GL_TRIANGLES);
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
void Renderer::renderMeshWithMaterial(RenderCall* rc)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material )
		return;
    assert(glGetError() == GL_NO_ERROR);

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	Camera* camera = Camera::current;
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

void Renderer::renderMeshWithMaterialLight(RenderCall* rc)
{
	//in case there is nothing to do
	if (!rc->mesh || !rc->mesh->getNumVertices() || !rc->material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (render_boundaries)
		rc->mesh->renderBounding(rc->model, true);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	Camera* camera = Camera::current;
	GFX::Texture* white = GFX::Texture::getWhiteTexture(); //a 1x1 white texture

	GFX::Texture* albedo_texture = rc->material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = rc->material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	GFX::Texture* metalic_texture = rc->material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	GFX::Texture* normal_map = rc->material->textures[SCN::eTextureChannel::NORMALMAP].texture;
	//GFX::Texture* occlusion_texture = rc->material->textures[SCN::eTextureChannel::OCCLUSION].texture;

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

	shader->setUniform("u_color", rc->material->color);
	shader->setUniform("u_emissive_factor", rc->material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white, 1);
	shader->setUniform("u_ambient", scene->ambient_light);
	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", rc->material->alpha_mode == SCN::eAlphaMode::MASK ? rc->material->alpha_cutoff : 0.001f);
	
	//set texture flags
	vec4 texture_flags = vec4(0,0,0,0); //normal, occlusion,specular | if 0 there is no texture 

	//normal map
	if (normal_map && enable_normal_map) {
		texture_flags.x = 1;
		shader->setUniform("u_normal_map",normal_map,2);
	}
	//occ & specular
	if (metalic_texture && (enable_occ|| enable_specular)) {
		if(enable_occ)
			texture_flags.y = 1;
		if (enable_specular) {
			texture_flags.z = 1;
			shader->setUniform("u_metalic_roughness", vec2(rc->material->metallic_factor, rc->material->roughness_factor));
			shader->setUniform("u_view_pos", Camera::current->eye);
			//std::cout << rc->material->metallic_factor << "\n";
		}

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

	glBlendFunc(GL_SRC0_ALPHA, GL_ONE);

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

#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ImGui::Checkbox("Wireframe", &render_wireframe);
	ImGui::Checkbox("Boundaries", &render_boundaries);

	//add here your stuff
	ImGui::Checkbox("Render Priority", &enable_render_priority);
	ImGui::Checkbox("Normal Map", &enable_normal_map);
	ImGui::Checkbox("Occlusion", &enable_occ);
	ImGui::Checkbox("Specular", &enable_specular);
	ImGui::Combo("Render Mode",(int*) &render_mode, "FLAT\0LIGHTS_MULTIPASS\0LIGHTS_SINGLEPASS\0", 3);
}

#else
void Renderer::showUI() {}
#endif