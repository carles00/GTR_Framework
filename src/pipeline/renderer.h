#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"

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
		LIGHTS
	};

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class RenderCall {
	public:
		GFX::Mesh* mesh;
		Material* material;
		Matrix44 model;
		
		static bool render_call_sort(RenderCall a, RenderCall b) {
			return (a.distance_to_camera > b.distance_to_camera );
		}
		//&& a.material->alpha_mode != SCN::eAlphaMode::BLEND && b.material->alpha_mode == SCN::eAlphaMode::BLEND

		float distance_to_camera;
	};

	class Renderer
	{
	public:
		bool render_wireframe;
		bool render_boundaries;
		eRenderMode render_mode;
		int N_LIGHTS;

		GFX::Texture* skybox_cubemap;

		SCN::Scene* scene;
		//render calls vector
		std::vector<RenderCall> render_order;
		std::vector<LightEntity*> lights;
		std::vector<LightEntity*> visible_lights;

		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene();

		//add here your functions
		//...
		void priorityRendering();
		void createRenderCall(Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);
		void walkEntities(SCN::Node* node, Camera* camera);
		void singlePass();
		void multiPass();

		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(RenderCall* rc);

		void renderMeshWithMaterialLight(RenderCall* rc);

		void showUI();

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader
	};

};