//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs
light_multipass basic.vs light_multipass.fs
light_singlepass basic.vs light_singlepass.fs
gbuffers basic.vs gbuffers.fs

deferred_global quad.vs deferred_global.fs
deferred_light quad.vs deferred_light.fs

\basic.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;
in vec4 a_color;

uniform vec3 u_camera_pos;

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

// Diffuse Reflections: Disney BRDF using retro-reflections using F term, this is much more complex!!
float Fd_Burley ( const in float NoV, const in float NoL,
const in float LoH, 
const in float linearRoughness)
{
        float f90 = 0.5 + 2.0 * linearRoughness * LoH * LoH;
		float lightScatter = F_Schlick(NoL, 1.0);
        float viewScatter  = F_Schlick(NoV, 1.0);      
        return lightScatter * viewScatter * RECIPROCAL_PI;
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
	
	vec3 f0 = mix(vec3(0.5f), albedo.xyz, metalic_roughness.y);
	vec3 diffuseColor = (1.0 - metalic_roughness.y) * albedo.xyz;

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
		vec3 Fr_d = specularBRDF(metalic_roughness.z, f0, NdotH, NdotV, NdotL, LdotH);
			
		float linearRoughness = metalic_roughness.z *metalic_roughness.z;
		vec3 Fd_d = diffuseColor * Fd_Burley(NdotV,NdotL,LdotH,linearRoughness); 
		
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
		light += max(NdotL, 0.0)* u_light_color * att_factor * shadow_factor * direct;		
	}
	else if(int(u_light_info.x) == DIRECTIONAL_LIGHT){
		float NdotL = dot(N, u_light_front);
		light += max(NdotL, 0.0)* u_light_color * shadow_factor;		
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
uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_emissive_factor;


#include "normal_functions"

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
	if(u_texture_flags.x == 1){
		vec3 normal_pixel = texture( u_normal_texture, v_uv ).xyz;
		N = perturbNormal(N,v_world_position, v_uv, normal_pixel);
	}
	

	vec3 emissive = u_emissive_factor * texture(u_emissive_texture, v_uv).xyz;
	vec3 metallicRoughness = texture(u_metalic_roughness, v_uv).xyz;

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
	float depth = texture( u_depth_texture, v_uv ).x;
	if(depth == 1.0)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_position = proj_worldpos.xyz / proj_worldpos.w;
	
	//gbuffers
	vec4 albedo = texture( u_albedo_texture, v_uv );
	vec4 extra = texture( u_extra_texture, v_uv );
	vec4 normal_info = texture( u_normal_texture, v_uv );
	vec4 metalic_roughness = texture(u_metalic_roughness, v_uv);

	vec3 N = normalize( normal_info.xyz * 2.0 - vec3(1.0) );
	vec3 V = normalize(u_eye - world_position);

	float shadow_factor = 1.0;
	if(u_shadow_params.x != 0.0){
		shadow_factor = testShadow(world_position);
	}

	//store light
	vec3 light = vec3(0.0);
	//we compute the reflection in base to the color and the metalness
	vec3 f0 = mix( vec3(0.5), albedo.xyz, metalic_roughness.x );

	//metallic materials do not have diffuse
	vec3 diffuseColor = (1.0 - metalic_roughness.x) * albedo.xyz;
	
	
	if(int(u_light_info.x) == POINT_LIGHT || int(u_light_info.x) == SPOT_LIGHT){
		vec3 L = u_light_position - world_position;
		float dist = length(L);
		L /= dist;
		float NdotL = dot(N, L);
		vec3 H = (V + L) / 2;
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);
		//compute the specular 
		vec3 Fr_d = specularBRDF(  metalic_roughness.y, f0, NdotH, NdotV, NdotL, LdotH);

		// Here we use the Burley, but you can replace it by the Lambert.
		float linearRoughness = metalic_roughness.y * metalic_roughness.y;
		vec3 Fd_d = diffuseColor * Fd_Burley(NdotV,NdotL,LdotH,linearRoughness); 

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
		float NdotL = dot(N, u_light_front);
		light += max(NdotL, 0.0)* u_light_color * shadow_factor;		
	}

	if(int(u_light_info.x) == NO_LIGHT){
		light = vec3(1.0);
	}

	vec4 color = vec4(0.0,0.0,0.0,1.0);
	color.xyz = light * albedo.xyz;

	FragColor = color;
	gl_FragDepth = depth;
}