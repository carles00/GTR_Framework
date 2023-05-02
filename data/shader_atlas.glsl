//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs
light_multipass basic.vs light_multipass.fs
light_singlepass basic.vs light_singlepass.fs

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

\light_multipass.fs

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


uniform vec4 u_light_info; //type, near, far, 0
uniform vec3 u_light_front;
uniform vec3 u_light_position;
uniform vec3 u_light_color;
uniform vec2 u_light_cone; // cos(max_ang), cos(min_ang)
uniform float u_max_distance;


out vec4 FragColor;

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
	
	
	if(int(u_light_info.x) == POINT_LIGHT || int(u_light_info.x) == SPOT_LIGHT){
		vec3 L = u_light_position - v_world_position;
		float dist = length(L);
		L /= dist;
		float NdotL = dot(N, L);

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
		light += max(NdotL, 0.0)* u_light_color * att_factor;		
	}
	else if(int(u_light_info.x) == DIRECTIONAL_LIGHT){
		float NdotL = dot(N, u_light_front);
		light += max(NdotL, 0.0)* u_light_color;

		if(u_texture_flags.z == 1){
			vec3 L = normalize(u_light_front);

			vec2 spec_factors = texture( u_occ_met_rough_texture, uv ).yz;
			vec3 V = normalize(u_view_pos - v_world_position);
			vec3 R = reflect(-L, N);
			float spec = pow(max(dot(V,R),0.0), max(u_metalic_roughness.x - spec_factors.x,0) );
			vec3 specular = (spec_factors.y - u_metalic_roughness.y)* spec * u_light_color;
			light += specular; 
		}
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

uniform vec4 u_light_info[MAX_LIGHTS]; //type, near, far, 0
uniform vec3 u_light_front[MAX_LIGHTS];
uniform vec3 u_light_position[MAX_LIGHTS];
uniform vec3 u_light_color[MAX_LIGHTS];
uniform vec2 u_light_cone[MAX_LIGHTS]; // cos(max_ang), cos(min_ang)
uniform int u_num_lights;

out vec4 FragColor;

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
				light += max(NdotL, 0.0)* u_light_color[i] * att_factor;		
			}
			else if(int(u_light_info[i].x) == DIRECTIONAL_LIGHT){
				float NdotL = dot(N, u_light_front[i]);
				light += max(NdotL, 0.0)* u_light_color[i];

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