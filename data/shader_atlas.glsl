//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs
light_multipass basic.vs light_multipass.fs
light_singlepass basic.vs light_singlepass.fs
gbuffers basic.vs gbuffers.fs
reflectionProbe basic.vs reflectionProbe.fs

tonemapper quad.vs tonemapper.fs
ssao quad.vs ssao.fs

spherical_probe basic.vs spherical_probe.fs

irradiance quad.vs irradiance.fs

deferred_global quad.vs deferred_global.fs
deferred_light quad.vs deferred_light.fs
deferred_ws basic.vs deferred_light.fs

volumetric quad.vs volumetric.fs

decal basic.vs decal.fs

fx_color quad.vs color_correction.fs
fx_blur quad.vs blur.fs
fx_motion_blur quad.vs motion_blur.fs

\basic.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;
in vec4 a_color;

uniform vec3 u_camera_position;

uniform mat4 u_model;
uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;
out vec4 v_color;

uniform float u_time;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( v_position, 1.0) ).xyz;
	
	//store the color in the varying var to use it from the pixel shader
	v_color = a_color;

	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\quad.vs

#version 330 core

in vec3 a_vertex;
in vec2 a_coord;
out vec2 v_uv;

void main()
{	
	v_uv = a_coord;
	gl_Position = vec4( a_vertex, 1.0 );
}


\flat.fs

#version 330 core

uniform vec4 u_color;

out vec4 FragColor;

void main()
{
	FragColor = u_color;
}


\texture.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );

	if(color.a < u_alpha_cutoff)
		discard;

	FragColor = color;
}


\skybox.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;

uniform samplerCube u_texture;
uniform vec3 u_camera_position;
out vec4 FragColor;

void main()
{
	vec3 E = v_world_position - u_camera_position;
	vec4 color = texture( u_texture, E );
	FragColor = color;
}


\multi.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	FragColor = color;
	NormalColor = vec4(N,1.0);
}


\depth.fs

#version 330 core

uniform vec2 u_camera_nearfar;
uniform sampler2D u_texture; //depth map
in vec2 v_uv;
out vec4 FragColor;

void main()
{
	float n = u_camera_nearfar.x;
	float f = u_camera_nearfar.y;
	float z = texture2D(u_texture,v_uv).x;
	if( n == 0.0 && f == 1.0 )
		FragColor = vec4(z);
	else
		FragColor = vec4( n * (z + 1.0) / (f + n - z * (f - n)) );
}


\instanced.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;

in mat4 u_model;

uniform vec3 u_camera_pos;

uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( a_vertex, 1.0) ).xyz;
	
	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\lights
#define NO_LIGHT 0
#define POINT_LIGHT 1
#define SPOT_LIGHT 2
#define DIRECTIONAL_LIGHT 3

uniform vec4 u_light_info; //type, near, far, 0
uniform vec3 u_light_front;
uniform vec3 u_light_position;
uniform vec3 u_light_color;
uniform vec2 u_light_cone; // cos(max_ang), cos(min_ang)

uniform vec2 u_shadow_params; //0 or 1, bias
uniform vec4 u_shadow_region;
uniform mat4 u_shadow_viewproj;
uniform sampler2D u_shadowmap;
float testShadow(vec3 pos)
{
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj * vec4(pos,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = proj_pos.xy / proj_pos.w;

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_params.y) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;

	//read depth from depth buffer in [0..+1] non-linear
	//float shadow_depth = texture( u_shadowmap, vec2(shadow_uv.x*0.25, shadow_uv.y*0.5)).x;
	float shadow_depth = texture( u_shadowmap, vec2((shadow_uv.x*u_shadow_region.z) + u_shadow_region.x, (shadow_uv.y*u_shadow_region.w)+u_shadow_region.y)).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if( shadow_depth < real_depth )
		shadow_factor = 0.0;
	//it is outside on the sides
	if( shadow_uv.x < 0.0 || shadow_uv.x > 1.0 ||
	    shadow_uv.y < 0.0 || shadow_uv.y > 1.0 )
			return 1.0;

	//it is before near or behind far plane
	if(real_depth < 0.0 || real_depth > 1.0)
		return 1.0;

	return shadow_factor;
}


\normal_functions

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
  // get edge vectors of the pixel triangle
  vec3 dp1 = dFdx(p);
  vec3 dp2 = dFdy(p);
  vec2 duv1 = dFdx(uv);
  vec2 duv2 = dFdy(uv);

  // solve the linear system
  vec3 dp2perp = cross(dp2, N);
  vec3 dp1perp = cross(N, dp1);
  vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
  vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

  // construct a scale-invariant frame 
  float invmax = 1.0 / sqrt(max(dot(T,T), dot(B,B)));
  return mat3(normalize(T * invmax), normalize(B * invmax), N);
}

vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255./127. - 128./127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}


\pbr_utils

#define RECIPROCAL_PI 0.3183098861837697
#define PI 3.14159265359

vec3 Fd_Lambert(vec3 albedo) {
    return albedo/PI;
}

// Fresnel term with scalar optimization(f90=1)
float F_Schlick( const in float VoH, 
const in float f0)
{
	float f = pow(1.0 - VoH, 5.0);
	return f0 + (1.0 - f0) * f;
}



// Fresnel term with colorized fresnel
vec3 F_Schlick( const in float VoH, 
const in vec3 f0)
{
	float f = pow(1.0 - VoH, 5.0);
	return f0 + (vec3(1.0) - f0) * f;
}


// Geometry Term: Geometry masking/shadowing due to microfacets
float GGX(float NdotV, float k){
	return NdotV / (NdotV * (1.0 - k) + k);
}

// Normal Distribution Function using GGX Distribution
float D_GGX (	const in float NoH, 
const in float linearRoughness )
{
	float a2 = linearRoughness * linearRoughness;
	float f = (NoH * NoH) * (a2 - 1.0) + 1.0;
	return a2 / (PI * f * f);
}

float G_Smith( float NdotV, float NdotL, float roughness)
{
	float k = pow(roughness + 1.0, 2.0) / 8.0;
	return GGX(NdotL, k) * GGX(NdotV, k);
}

//this is the cook torrance specular reflection model
vec3 specularBRDF( float roughness, vec3 f0, 
float NoH, float NoV, float NoL, float LoH )
{
	float a = roughness * roughness;

	// Normal Distribution Function
	float D = D_GGX( NoH, a );

	// Fresnel Function
	vec3 F = F_Schlick( LoH, f0 );

	// Visibility Function (shadowing/masking)
	float G = G_Smith( NoV, NoL, roughness );
		
	// Norm factor
	vec3 spec = D * G * F;
	spec /= (4.0 * NoL * NoV + 1e-6);

	return spec;
}


\light_multipass.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform vec3 u_emissive_factor;
uniform vec2 u_metalic_roughness; //metalic, roughness
uniform vec3 u_view_pos;

uniform vec4 u_texture_flags; //normal, occlusion, specular, pbr
uniform sampler2D u_albedo_texture;
uniform sampler2D u_emissive_texture;
uniform sampler2D u_occ_met_rough_texture;
uniform sampler2D u_normal_map;

uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_ambient;

#include "lights"
#include "normal_functions"
#include "pbr_utils"

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 albedo = u_color;
	float shadow_factor = 1.0;
	if(u_shadow_params.x != 0.0){
		shadow_factor = testShadow(v_world_position);
	}

	albedo *= texture( u_albedo_texture, v_uv );
	//discard if alpha is too low
	if(albedo.a < u_alpha_cutoff)
		discard;

	//normal
	vec4 metalic_roughness = texture(u_occ_met_rough_texture, uv);
	metalic_roughness.g = pow(metalic_roughness.g, u_metalic_roughness.x);
	metalic_roughness.b = pow(metalic_roughness.b, u_metalic_roughness.y);
	vec3 N = normalize(v_normal);
	vec3 V = normalize(u_view_pos - v_world_position);

	//if the mesh has normal map
	if(u_texture_flags.x == 1){ 
		vec3 normal_pixel = texture( u_normal_map, uv ).xyz;
		N = perturbNormal(N,v_world_position,uv , normal_pixel);
	}

	
	//store light
	vec3 light = vec3(0.0);
	//add ambient
	//occulision enabled
	if(u_texture_flags.y == 1){
		float occ_fact = metalic_roughness.x;
		light += u_ambient * occ_fact;
	}else{
		light += u_ambient;
	}
	
	

	if(int(u_light_info.x) == POINT_LIGHT || int(u_light_info.x) == SPOT_LIGHT){
		vec3 L = u_light_position - v_world_position;
		float dist = length(L);
		L /= dist;
		vec3 H = (V + L) / 2;
		float NdotV = dot(N, V);
		float NdotL = dot(N, L);
		float NdotH = dot(N, H);
		float LdotH = dot(L, H);

		//pbr
		vec3 f0 = mix(vec3(0.5f), albedo.xyz, metalic_roughness.b);
		vec3 diffuseColor = (1.0 - metalic_roughness.b) * albedo.xyz;
		vec3 Fr_d = specularBRDF(metalic_roughness.g, f0, NdotH, NdotV, NdotL, LdotH);
			
		float linearRoughness = metalic_roughness.b *metalic_roughness.b;
		vec3 Fd_d = diffuseColor * Fd_Lambert(albedo.xyz); 
		
		vec3 direct = Fr_d + Fd_d;

		//attenuation
		float att_factor = (u_light_info.z - dist) / u_light_info.z;
		att_factor = max(att_factor, 0);

		if(int(u_light_info.x) == SPOT_LIGHT){
			float cos_angle = dot(u_light_front, L);
			if(cos_angle < u_light_cone.y){
				att_factor = 0;
			}
			else if(cos_angle < u_light_cone.x){
				att_factor *= 1.0 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
			}
		}
		if(u_texture_flags.w == 0)
			direct = vec3(1.0,1.0,1.0);
		vec3 light_params = max(NdotL, 0.0)* u_light_color * att_factor * shadow_factor;
		light += light_params * direct;
	}
	else if(int(u_light_info.x) == DIRECTIONAL_LIGHT){
		vec3 L = u_light_front;
		float NdotL = dot(N, L);
		vec3 H = (V + L) / 2;
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);
		vec3 f0 = mix(vec3(0.5f), albedo.xyz, metalic_roughness.b);
		vec3 diffuseColor = (1.0 - metalic_roughness.b) * albedo.xyz;
		vec3 Fr_d = specularBRDF(metalic_roughness.g, f0, NdotH, NdotV, NdotL, LdotH);
			
		float linearRoughness = metalic_roughness.b *metalic_roughness.b;
		vec3 Fd_d = diffuseColor * Fd_Lambert(albedo.xyz); 
		
		vec3 direct = Fr_d + Fd_d;



		light += max(NdotL, 0.0)* u_light_color * shadow_factor * direct;		
	}
	
	//specular
	if(u_texture_flags.z == 1){ 
		vec3 L = vec3(0,0,0);
		if(int(u_light_info.x) == DIRECTIONAL_LIGHT){
			L = normalize(u_light_front);
		}else{
			L = normalize( u_light_position - v_world_position);
		}

		vec2 spec_factors = texture( u_occ_met_rough_texture, uv ).yz;
		vec3 R = reflect(-L, N);
		float spec = pow(max(dot(V,R),0.0), max(u_metalic_roughness.x - spec_factors.x,0) );
		vec3 specular = (spec_factors.y - u_metalic_roughness.y)* spec * u_light_color;
		light += specular; 
	}
	if(int(u_light_info.x) == NO_LIGHT){
		light = vec3(1.0);
	}

	vec3 color = albedo.xyz * light;
	color += u_emissive_factor * texture( u_emissive_texture, v_uv ).xyz;

	FragColor = vec4(color, albedo.a);
}

\light_singlepass.fs

#version 330 core

#define NO_LIGHT 0
#define POINT_LIGHT 1
#define SPOT_LIGHT 2
#define DIRECTIONAL_LIGHT 3



in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform vec3 u_emissive_factor;
uniform vec2 u_metalic_roughness; //metalic, roughness
uniform vec3 u_view_pos;

uniform vec4 u_texture_flags; //normal, occlusion, specular
uniform sampler2D u_albedo_texture;
uniform sampler2D u_emissive_texture;
uniform sampler2D u_occ_met_rough_texture;
uniform sampler2D u_normal_map;

uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_ambient;

const int MAX_LIGHTS = 10;
//lights
uniform vec4 u_light_info[MAX_LIGHTS]; //type, near, far, 0
uniform vec3 u_light_front[MAX_LIGHTS];
uniform vec3 u_light_position[MAX_LIGHTS];
uniform vec3 u_light_color[MAX_LIGHTS];
uniform vec2 u_light_cone[MAX_LIGHTS]; // cos(max_ang), cos(min_ang)
uniform int u_num_lights;

//shadows
uniform vec2 u_shadow_params[MAX_LIGHTS];
uniform vec4 u_shadow_region[MAX_LIGHTS];
uniform mat4 u_shadow_viewproj[MAX_LIGHTS];
uniform sampler2D u_shadowmap;

#include "normal_functions"

out vec4 FragColor;

float testShadow(vec3 pos, int i)
{
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj[i] * vec4(pos,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = proj_pos.xy / proj_pos.w;

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_params[i].y) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;

	//read depth from depth buffer in [0..+1] non-linear
	//float shadow_depth = texture( u_shadowmap, vec2(shadow_uv.x*0.25, shadow_uv.y*0.5)).x;
	float shadow_depth = texture( u_shadowmap, vec2((shadow_uv.x*u_shadow_region[i].z) + u_shadow_region[i].x, (shadow_uv.y*u_shadow_region[i].w)+u_shadow_region[i].y)).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if( shadow_depth < real_depth )
		shadow_factor = 0.0;
	//it is outside on the sides
	if( shadow_uv.x < 0.0 || shadow_uv.x > 1.0 ||
	    shadow_uv.y < 0.0 || shadow_uv.y > 1.0 )
			return 1.0;

	//it is before near or behind far plane
	if(real_depth < 0.0 || real_depth > 1.0)
		return 1.0;

	return shadow_factor;
}


void main()
{
	vec2 uv = v_uv;
	vec4 albedo = u_color;
	albedo *= texture( u_albedo_texture, v_uv );
	//discard if alpha is too low
	if(albedo.a < u_alpha_cutoff)
		discard;

	//normal
	vec3 N = normalize(v_normal);
	
	//if the mesh has normal map
	if(u_texture_flags.x == 1){ 
		vec3 normal_pixel = texture( u_normal_map, uv ).xyz;
		N = perturbNormal(N,v_world_position,uv , normal_pixel);
	}

	
	//store light
	vec3 light = vec3(0.0);
	//add ambient
	//occulision enabled
	if(u_texture_flags.y == 1){
		float occ_fact = texture( u_occ_met_rough_texture, uv ).x;
		light += u_ambient * occ_fact;
	}else{
		light += u_ambient;
	}
	
	for(int i = 0; i<MAX_LIGHTS; i++)
	{
		if(i < u_num_lights){
			float shadow_factor = 1.0;
			if(u_shadow_params[i].x != 0.0){
				shadow_factor = testShadow(v_world_position,i);
			}

			if(int(u_light_info[i].x) == POINT_LIGHT || int(u_light_info[i].x) == SPOT_LIGHT){
				vec3 L = u_light_position[i] - v_world_position;
				float dist = length(L);
				L /= dist;
				float NdotL = dot(N, L);

				//attenuation
				float att_factor = (u_light_info[i].z - dist) / u_light_info[i].z;
				att_factor = max(att_factor, 0);

				if(int(u_light_info[i].x) == SPOT_LIGHT){
					float cos_angle = dot(u_light_front[i], L);
					if(cos_angle < u_light_cone[i].y){
						att_factor = 0;
					}
					else if(cos_angle < u_light_cone[i].x){
						att_factor *= 1.0 - (cos_angle - u_light_cone[i].x) / (u_light_cone[i].y - u_light_cone[i].x);
					}
				}
				light += max(NdotL, 0.0)* u_light_color[i] * att_factor* shadow_factor;		
			}
			else if(int(u_light_info[i].x) == DIRECTIONAL_LIGHT){
				float NdotL = dot(N, u_light_front[i]);
				light += max(NdotL, 0.0)* u_light_color[i] * shadow_factor;

				if(u_texture_flags.z == 1){
					vec3 L = normalize(u_light_front[i]);

					vec2 spec_factors = texture( u_occ_met_rough_texture, uv ).yz;
					vec3 V = normalize(u_view_pos - v_world_position);
					vec3 R = reflect(-L, N);
					float spec = pow(max(dot(V,R),0.0), max(u_metalic_roughness.x - spec_factors.x,0) );
					vec3 specular = (spec_factors.y - u_metalic_roughness.y)* spec * u_light_color[i];
					light += specular; 
				}
			}	
		}
	}

	vec3 color = albedo.xyz * light;
	color += u_emissive_factor * texture( u_emissive_texture, v_uv ).xyz;

	FragColor = vec4(color, albedo.a);
}

\gbuffers.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform vec4 u_texture_flags; //normal, occlusion, specular
uniform sampler2D u_albedo_texture;
uniform sampler2D u_emissive_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_metalic_roughness;
uniform vec3 u_camera_position;

uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_emissive_factor;
uniform vec2 u_metalness_roughness;

#include "normal_functions"

uniform samplerCube u_environment;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;
layout(location = 2) out vec4 ExtraColor;
layout(location = 3) out vec4 MetalRoughColor;


void main()
{
	vec4 color = u_color;
	color *= texture( u_albedo_texture, v_uv );
	
	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);
	vec3 E = normalize(v_world_position - u_camera_position);
	vec3 R = reflect(E, N);

	if(u_texture_flags.x == 1){
		vec3 normal_pixel = texture( u_normal_texture, v_uv ).xyz;
		N = perturbNormal(N,v_world_position, v_uv, normal_pixel);
	}
	//occlusion
	float occ_fact = 1.0;
	if(u_texture_flags.y == 1){
		float occ_fact = texture( u_metalic_roughness, v_uv ).x;
	}

	vec3 emissive = u_emissive_factor * texture(u_emissive_texture, v_uv).xyz;
	vec3 metallicRoughness = texture(u_metalic_roughness, v_uv).xyz;


	//float fresnel = 1.0 - max(0.0, dot(-E, N));
	float reflective_factor =  metallicRoughness.b * u_metalness_roughness.y;

	vec3 reflected_color = texture(u_environment, R).xyz;
	color.xyz = mix(color.xyz, reflected_color, reflective_factor);
	metallicRoughness.g = pow(metallicRoughness.g, u_metalness_roughness.x);
	metallicRoughness.b = pow(metallicRoughness.b, u_metalness_roughness.y);

	FragColor = vec4(color.xyz, 1.0);
	NormalColor = vec4(N*0.5 + vec3(0.5),1.0);
	ExtraColor = vec4(emissive, 1.0);
	MetalRoughColor = vec4(metallicRoughness,1.0);

}

\deferred_global.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_extra_texture;
uniform sampler2D u_depth_texture;

uniform vec3 u_ambient_light;

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;

	float depth = texture( u_depth_texture, v_uv ).x;
	if(depth == 1.0)
		discard;

	vec4 albedo = texture( u_albedo_texture, v_uv );
	vec4 extra = texture( u_extra_texture, v_uv );
	//vec4 normal_info = texture( u_normal_texture, v_uv );
	//vec3 N = normalize( normal_info.xyz * 2.0 - vec3(1.0) );

	vec4 color = vec4(0.0);

	color.xyz += extra.xyz + u_ambient_light * albedo.xyz;
	FragColor = color;
	gl_FragDepth = depth;
}


\deferred_light.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_extra_texture;
uniform sampler2D u_depth_texture;
uniform sampler2D u_metalic_roughness;

uniform mat4 u_ivp;
uniform vec2 u_iRes;
uniform vec3 u_eye;
uniform float u_pbr_state;


#include "lights"
#include "pbr_utils"

out vec4 FragColor;

void main()
{
	
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;
	if(depth == 1.0)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;
	
	//gbuffers
	vec4 albedo = texture( u_albedo_texture, uv );
	vec4 extra = texture( u_extra_texture, uv );
	vec4 normal_info = texture( u_normal_texture, uv );
	vec4 metalic_roughness = texture(u_metalic_roughness, uv);

	vec3 N = normalize( normal_info.xyz * 2.0 - vec3(1.0) );
	vec3 V = normalize(u_eye - world_position);

	float shadow_factor = 1.0;
	if(u_shadow_params.x != 0.0){
		shadow_factor = testShadow(world_position);
	}

	//store light
	vec3 light = vec3(0.0);
	
	
	if(int(u_light_info.x) == POINT_LIGHT || int(u_light_info.x) == SPOT_LIGHT){
		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix( vec3(0.5), albedo.xyz, metalic_roughness.b );

		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalic_roughness.b) * albedo.xyz;
		vec3 L = u_light_position - world_position;
		float dist = length(L);
		L /= dist;
		float NdotL = dot(N, L);
		vec3 H = (V + L) / 2;
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);
		//compute the specular 
		vec3 Fr_d = specularBRDF(  metalic_roughness.g, f0, NdotH, NdotV, NdotL, LdotH);

		// Here we use the Burley, but you can replace it by the Lambert.
		vec3 Fd_d = diffuseColor * Fd_Lambert(albedo.xyz); 

		//add diffuse and specular reflection
		vec3 direct = Fr_d + Fd_d;
		//attenuation
		float att_factor = (u_light_info.z - dist) / u_light_info.z;
		att_factor = max(att_factor, 0);

		if(int(u_light_info.x) == SPOT_LIGHT){
			float cos_angle = dot(u_light_front, L);
			if(cos_angle < u_light_cone.y){
				att_factor = 0;
			}
			else if(cos_angle < u_light_cone.x){
				att_factor *= 1.0 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
			}
		}
		if(u_pbr_state == 0.0)
			direct = vec3(1.0,1.0,1.0);
		light += max(NdotL, 0.0)* u_light_color * att_factor * shadow_factor * direct;
		 
				
	}
	else if(int(u_light_info.x) == DIRECTIONAL_LIGHT){
		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix( vec3(0.5), albedo.xyz, metalic_roughness.b );

		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalic_roughness.b) * albedo.xyz; 
		vec3 L = u_light_front;
		float NdotL = dot(N, L);
		vec3 H = (V + L) / 2;
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);
		//compute the specular 
		vec3 Fr_d = specularBRDF(  metalic_roughness.g, f0, NdotH, NdotV, NdotL, LdotH);

		// Here we use the Burley, but you can replace it by the Lambert.
		vec3 Fd_d = diffuseColor * Fd_Lambert(albedo.xyz); 

		//add diffuse and specular reflection
		vec3 direct = Fr_d + Fd_d;

		if(u_pbr_state == 0.0)
			direct = vec3(1.0,1.0,1.0);
			
		light += max(NdotL, 0.0)* u_light_color * shadow_factor * direct;		
	}

	if(int(u_light_info.x) == NO_LIGHT){
		light = vec3(1.0);
	}

	vec4 color = vec4(0.0,0.0,0.0,1.0);
	color.xyz = light * albedo.xyz *extra.w;

	FragColor = color;
	gl_FragDepth = depth;
}

\tonemapper.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_texture;

uniform float u_scale; //color scale before tonemapper
uniform float u_average_lum; 
uniform float u_lumwhite2;
uniform float u_igamma; //inverse gamma

out vec4 FragColor;

void main() {
	vec4 color = texture2D( u_texture, v_uv );
	vec3 rgb = color.xyz;

	float lum = dot(rgb, vec3(0.2126, 0.7152, 0.0722));
	float L = (u_scale / u_average_lum) * lum;
	float Ld = (L * (1.0 + L / u_lumwhite2)) / (1.0 + L);

	rgb = (rgb / lum) * Ld;
	rgb = max(rgb,vec3(0.001));
	rgb = pow( rgb, vec3( u_igamma ) );
	FragColor = vec4( rgb, color.a );
}




\ssao.fs

#version 330 core

#define NUM_POINTS 64

in vec2 v_uv;

uniform sampler2D u_normal_texture;
uniform sampler2D u_depth_texture;


uniform mat4 u_viewprojection;
uniform mat4 u_ivp;
uniform vec2 u_iRes;
uniform vec3 u_random_points[NUM_POINTS];
uniform float u_radius;
uniform float u_ssao_plus;


#include "normal_functions"

layout(location = 0) out vec4 FragColor;


void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;
	vec4 normal_text = texture( u_normal_texture, uv);
	vec3 N = normalize( normal_text.xyz * 2.0 - vec3(1.0) );
	float ao = 1.0;
	if(depth < 1.0)
	{
		vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
		vec4 proj_worldpos = u_ivp * screen_pos;
		vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;
		
		//lets use 64 samples
		const int samples = NUM_POINTS;
		int num = samples; //num samples that passed the are outside

		//for every sample around the point
		for( int i = 0; i < samples; ++i )
		{
			//compute is world position using the random
			vec3 p = world_position + u_random_points[i] * u_radius;
			
			if(u_ssao_plus == 1.0)
			{
				vec3 pointVector = u_random_points[i];
				if( dot(N, pointVector) < 0)
					p -= pointVector;
			}
			//find the uv in the depth buffer of this point
			vec4 proj = u_viewprojection * vec4(p,1.0);
			proj.xy /= proj.w; //convert to clipspace from homogeneous
			//apply a tiny bias to its z before converting to clip-space
			proj.z = (proj.z - 0.005) / proj.w;
			proj.xyz = proj.xyz * 0.5 + vec3(0.5); //to [0..1]
			//read p true depth
			float pdepth = texture( u_depth_texture, proj.xy ).x;
			//compare true depth with its depth
			if( pdepth < proj.z ) //if true depth smaller, is inside
				num--; //remove this point from the list of visible
		}
		ao = float(num) / float(samples);

	}
	ao = pow(ao, 1.0/2.2);
	FragColor = vec4(ao,ao,ao,1.0);
}

\probes_utils

const float Pi = 3.141592654;
const float CosineA0 = Pi;
const float CosineA1 = (2.0 * Pi) / 3.0;
const float CosineA2 = Pi * 0.25;
struct SH9 { float c[9]; }; //to store weights
struct SH9Color { vec3 c[9]; }; //to store colors

void SHCosineLobe(in vec3 dir, out SH9 sh) //SH9
{
	// Band 0
	sh.c[0] = 0.282095 * CosineA0;
	// Band 1
	sh.c[1] = 0.488603 * dir.y * CosineA1; 
	sh.c[2] = 0.488603 * dir.z * CosineA1;
	sh.c[3] = 0.488603 * dir.x * CosineA1;
	// Band 2
	sh.c[4] = 1.092548 * dir.x * dir.y * CosineA2;
	sh.c[5] = 1.092548 * dir.y * dir.z * CosineA2;
	sh.c[6] = 0.315392 * (3.0 * dir.z * dir.z - 1.0) * CosineA2;
	sh.c[7] = 1.092548 * dir.x * dir.z * CosineA2;
	sh.c[8] = 0.546274 * (dir.x * dir.x - dir.y * dir.y) * CosineA2;
}

vec3 ComputeSHIrradiance(in vec3 normal, in SH9Color sh)
{
	// Compute the cosine lobe in SH, oriented about the normal direction
	SH9 shCosine;
	SHCosineLobe(normal, shCosine);
	// Compute the SH dot product to get irradiance
	vec3 irradiance = vec3(0.0);
	for(int i = 0; i < 9; ++i)
		irradiance += sh.c[i] * shCosine.c[i];

	return irradiance;
}


\spherical_probe.fs

#version 330 core

in vec3 v_world_position;
in vec3 v_normal;

#include probes_utils
uniform vec3 u_coeffs[9];

out vec4 FragColor;

void main()
{
	vec4 color = vec4(1.0);
	vec3 N = normalize(v_normal);

	SH9Color sh;
	for (int i = 0; i < 9; ++i)
		sh.c[i] = u_coeffs[i];

	color.xyz = max(vec3(0.0), ComputeSHIrradiance(N, sh));

	FragColor = color;
}


\irradiance.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_depth_texture;

uniform vec3 u_irr_start;
uniform vec3 u_irr_end;
uniform vec3 u_irr_dims;
uniform float u_irr_normal_distance;
uniform float u_irr_multiplier;
uniform vec3 u_irr_delta;
uniform int u_num_probes;
uniform sampler2D u_probes_texture;

uniform mat4 u_ivp;
uniform vec2 u_iRes;

#include probes_utils

out vec4 FragColor;

void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;
	if(depth == 1.0)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 worldpos = proj_worldpos.xyz / proj_worldpos.w;
	
	//gbuffers
	vec4 albedo = texture( u_albedo_texture, uv );
	vec4 normal_info = texture( u_normal_texture, uv );

	vec3 N = normalize( normal_info.xyz * 2.0 - vec3(1.0) );
	
	//computing nearest probe index based on world position
	vec3 irr_range = u_irr_end - u_irr_start;
	vec3 irr_local_pos = clamp( worldpos - u_irr_start 
	+ N * u_irr_normal_distance, //offset a little
	vec3(0.0), irr_range );

	//convert from world pos to grid pos
	vec3 irr_norm_pos = irr_local_pos / u_irr_delta;

	//round values as we cannot fetch between rows for now
	vec3 local_indices = round( irr_norm_pos );

	//compute in which row is the probe stored
	float row = local_indices.x + 
	local_indices.y * u_irr_dims.x + 
	local_indices.z * u_irr_dims.x * u_irr_dims.y;

	//find the UV.y coord of that row in the probes texture
	float row_uv = (row + 1.0) / (u_num_probes + 1.0);

	SH9Color sh;

	//fill the coefficients
	const float d_uvx = 1.0 / 9.0;
	for(int i = 0; i < 9; ++i)
	{
		vec2 coeffs_uv = vec2( (float(i)+0.5) * d_uvx, row_uv );
		sh.c[i] = texture( u_probes_texture, coeffs_uv).xyz;
	}

	//now we can use the coefficients to compute the irradiance
	vec3 irradiance = max(vec3(0.0), ComputeSHIrradiance( N, sh ) * u_irr_multiplier);



	FragColor = vec4(albedo.xyz * irradiance,1.0);
}


\volumetric.fs

#version 330 core

#define SAMPLES 64

in vec2 v_uv;

uniform sampler2D u_depth_texture;
uniform mat4 u_ivp;
uniform vec2 u_iRes;
uniform vec3 u_camera_position;
uniform float u_air_density;
uniform vec3 u_ambient;
uniform float u_random;
uniform float u_time;

#include "lights"
#include "normal_functions"

layout(location = 0) out vec4 FragColor;


vec3 computeLight(vec3 pos)
{
	vec3 light = vec3(0.0);
	float shadow_factor = 1.0;

	if(u_shadow_params.x != 0.0)
		shadow_factor = testShadow(pos);

	if(int(u_light_info.x) == DIRECTIONAL_LIGHT)
	{
		light = u_light_color * shadow_factor;
	}
	else if (int(u_light_info.x) == POINT_LIGHT || int(u_light_info.x) == SPOT_LIGHT)
	{
		vec3 L = u_light_position - pos;
		float dist = length(L);
		L /= dist;
		float att = max(0.0, (u_light_info.z - dist)/u_light_info.z);
		if(int(u_light_info.x) == SPOT_LIGHT){
			float cos_angle = dot(u_light_front, L);
			if(cos_angle < u_light_cone.y){
				att = 0;
			}
			else if(cos_angle < u_light_cone.x){
				att *= 1.0 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
			}
		}

		light += u_light_color * att * shadow_factor;
	}
	return light;
}

float rand(vec2 co)
{
	return fract(sin(dot(co, vec2(12.9898, 78.233)))*43758.5453123);
}

void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;	
		
	//compute ray info
	vec3 ray_start = u_camera_position;
	vec3 ray_dir = ( world_position - ray_start );
	float ray_length = length(ray_dir);
	ray_dir /= ray_length;
	ray_dir = normalize(ray_dir);
	ray_length = min( 500.0, ray_length ); //max ray
	float step_dist = ray_length / float(SAMPLES);

	ray_start += ray_dir*rand(uv + vec2(u_random, u_time)) * step_dist;

	vec3 current_pos = ray_start;
	vec3 ray_offset = ray_dir * step_dist;

	vec3 volumetric = vec3(0.0);

	vec3 irradiance = vec3(0.0);
	float transparency = 1.0;
	float air_step = u_air_density * step_dist;

	for(int i = 0; i < SAMPLES; ++i)
	{
		//evaluate contribution
		vec3 light = computeLight(current_pos);

		//accumulate the amount of light
		irradiance += (u_ambient + light )* transparency * air_step;

		//advance to next position
		current_pos.xyz += ray_offset;

		transparency -= air_step;

		//too dense, nothing can be seen behind
		if( transparency < 0.001 )
			break;
	}

	FragColor = vec4(irradiance, 1.0 - clamp(transparency, 0.0, 1.0));
}

\decal.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_depth_texture;
uniform sampler2D u_color_texture;
uniform mat4 u_ivp;
uniform mat4 u_imodel;
uniform vec2 u_iRes;

uniform vec3 u_ambient_light;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;
layout(location = 2) out vec4 ExtraColor;
layout(location = 3) out vec4 MetalRoughColor;

void main()
{
	
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;
	if(depth == 1.0)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;
	
	vec3 localpos = (u_imodel * vec4(world_position, 1.0)).xyz ;
	
	if( localpos.x < -0.5 || localpos.x > 0.5 ||
    		localpos.y < -0.5 || localpos.y > 0.5 ||
    		localpos.z < -0.5 || localpos.z > 0.5 )
			discard;

	vec2 decal_uv = localpos.xy + vec2(0.5);
	
	vec4 color = texture(u_color_texture, decal_uv);

	FragColor = color;
	NormalColor  = vec4(0.0);
	ExtraColor = vec4(0.0);
	MetalRoughColor = vec4(0.0);
}

\reflectionProbe.fs

#version 330 core

in vec3 v_position;
in vec3 v_normal;
in vec3 v_world_position;

uniform samplerCube u_texture;
uniform vec3 u_camera_position;
out vec4 FragColor;

void main()
{
	vec3 N = normalize(v_normal);
	vec3 E = v_world_position - u_camera_position;
	vec3 R = reflect(E,N);
	vec4 color = textureLod( u_texture, R, 2.0);
	FragColor = color;
}


\color_correction.fs

#version 330 core

uniform sampler2D u_texture;
uniform float u_brightness;
in vec2 v_uv;
out vec4 FragColor;



void main()
{
	vec4 color = texture(u_texture, v_uv);

	color.xyz *= u_brightness;

	FragColor = color;
}

\blur.fs

#version 330 core

//linear blur shader of 9 samples
precision highp float;
uniform sampler2D u_texture;
in vec2 v_uv;
out vec4 FragColor;

uniform vec2 u_offset;
uniform float u_intensity;

void main() {
   vec4 sum = vec4(0.0);
   sum += texture(u_texture, v_uv + u_offset * -4.0) * 0.05/0.98;
   sum += texture(u_texture, v_uv + u_offset * -3.0) * 0.09/0.98;
   sum += texture(u_texture, v_uv + u_offset * -2.0) * 0.12/0.98;
   sum += texture(u_texture, v_uv + u_offset * -1.0) * 0.15/0.98;
   sum += texture(u_texture, v_uv) * 0.16/0.98;
   sum += texture(u_texture, v_uv + u_offset * 4.0) * 0.05/0.98;
   sum += texture(u_texture, v_uv + u_offset * 3.0) * 0.09/0.98;
   sum += texture(u_texture, v_uv + u_offset * 2.0) * 0.12/0.98;
   sum += texture(u_texture, v_uv + u_offset * 1.0) * 0.15/0.98;
   FragColor = sum * u_intensity;
}


\motion_blur.fs

#version 330 core

uniform sampler2D u_texture;
uniform sampler2D u_depth_texture;
uniform mat4 u_ivp;
uniform mat4 u_prev_vp;
uniform vec2 u_iRes;
in vec2 v_uv;
out vec4 FragColor;



void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture( u_depth_texture, uv ).x;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;

	vec4 prev_screenpos = u_prev_vp * vec4( world_position, 1.0 );
	prev_screenpos.xyz /= prev_screenpos.w;
	vec2 prev_uv = prev_screenpos.xy * 0.5 + vec2(0.5);
	
	vec4 color = vec4(0.0);
	for(int i = 0; i < 16; ++i)
	{
		vec2 int_uv = mix(uv, prev_uv, float(i)/16.0);
		color += texture(u_texture, int_uv);
	}

	color /= 16.0;

	FragColor = color;
}