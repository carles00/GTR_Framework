#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"
#include "../gfx/sphericalharmonics.h"

#define MAX_LIGHTS 10
//forward declarations

class Camera;
class Skeleton;
namespace GFX {
	class Shader;
	class Mesh;
	class FBO;
}

namespace SCN {

	class Prefab;
	class Material;
	enum eRenderMode {
		FLAT,
		TEXTURED,
		LIGHTS_MULTIPASS,
		LIGHTS_SINGLEPASS,
		DEFERRED
	};
	
	struct sLightsContainer {
		std::vector<LightEntity*> lights;
	};


	struct sProbe {
		vec3 pos; //where is located
		vec3 local; //its ijk pos in the matrix
		int index; //its index in the linear array
		SphericalHarmonics sh; //coeffs
	};

	struct sIrradianceCacheInfo {
		int num_probes;
		vec3 dims;
		vec3 start;
		vec3 end;
	};

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class RenderCall {
	public:
		GFX::Mesh* mesh;
		Material* material;
		Matrix44 model;
		
		static bool render_call_sort(RenderCall a, RenderCall b) {
			return (a.distance_to_camera < b.distance_to_camera);
		}
		static bool render_call_alpha_sort(RenderCall a, RenderCall b) {
			return (a.distance_to_camera > b.distance_to_camera);
		}


		float distance_to_camera;
	};

	class Renderer
	{
	public:
		//ImGui options
		bool render_wireframe;
		bool render_boundaries;
		bool show_shadowmaps;
		eRenderMode render_mode;
		bool enable_normal_map;
		bool enable_occ;
		bool enable_specular;
		bool enable_shadows;
		bool show_gbuffers;
		bool pbr_is_active;
		bool show_ssao;
		bool show_only_fbo;
		bool ssao_plus;
		bool swap_ssao;
		bool show_probes;
		bool update_probes;
		int buffers_to_show[4];
		float ssao_radius;
		std::vector<vec3> random_points;
		std::vector<vec3> copy_random_points;

		//shadows
		GFX::FBO* shadow_atlas_fbo;
		GFX::Texture* shadow_atlas;
		GFX::Texture* ssao_blur;
		GFX::Texture* probes_texture;

		float shadow_atlas_width;
		float shadow_atlas_height;
		float shadowmap_width;
		float shadowmap_height;

		float tonemapper_scale;
		float average_lum;
		float lumwhite2;
		float gamma;
		float irradiance_multiplier;

		GFX::Texture* skybox_cubemap;
		//deferred
		GFX::FBO* gbuffers_fbo;
		GFX::FBO* illumination_fbo;
		GFX::FBO* ssao_fbo;
		GFX::FBO* irr_fbo;

		sIrradianceCacheInfo irradiance_cache_info;
		std::vector<sProbe> probes;


		SCN::Scene* scene;
		//render calls vector
		std::vector<RenderCall> render_order;
		std::vector<RenderCall> render_order_alpha;
		std::vector<LightEntity*> lights;

		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene(Camera* camera);

		//add here your functions
		//...
		void renderFrame(Camera* camera);
		void captureIrradiance();
		void loadIrradianceCache();
		void uploadIrradianceCache();
		void createRenderCall(Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, vec3 camera_pos);
		void renderDeferred(Camera* camera);
		void renderForward(Camera* camera);
		void walkEntities(SCN::Node* node, Camera* camera);
		void singlePass(RenderCall* rc, GFX::Shader* shader);
		void multiPass(RenderCall* rc, GFX::Shader* shader);
		void generateShadowmaps();
		bool cullLights(LightEntity* light, BoundingBox bb);
		bool spotLightAABB(LightEntity* light, BoundingBox bb);
		void gbuffersToShader(GFX::FBO* gbuffers, GFX::Shader* shader);


		void renderMeshWithMaterialGBuffers(RenderCall* rc, Camera* camera);
		void renderAlphaObjects(Camera* camera);

		//debug
		void showShadowmaps();
		
		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);

		void renderSceneNodes(Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera);

		void renderProbe(sProbe& probe);
		void applyIrradiance();
		void captureProbe(sProbe& probe);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(RenderCall* rc, Camera* camera);
		void renderMeshWithMaterialFlat(RenderCall* rc, Camera* camera);
		void renderMeshWithMaterialLight(RenderCall* rc, Camera* camera);

		void showUI();

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader
		void lightToShader(LightEntity* light, GFX::Shader* shader);
		std::vector<vec3> generateSpherePoints(int num, float radius, bool hemi);
	};

};

