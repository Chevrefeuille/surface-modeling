// Version d'OpenGL
#version 150

in vec3 vert_color;

in vec3 vert_normal;
in vec3 light_dir;

out vec4 frag_color;

//uniform sampler2D texSampler;

// Fonction appellee pour chaque fragment
void main()
{
  vec3 normal = normalize(vert_normal);

  float t = dot(-light_dir, normal);
/*  t = max(t, 0.0);*/

  vec3 color = vert_color;
  if(t < 0)
  {
    color = color * 0.1;
    t = -t;
  }
  color = color * (t*0.9 + 0.1);
  
  frag_color = vec4(color, 1.0);
//  frag_color = vec4(vert_normal, 1.0);
  
}
