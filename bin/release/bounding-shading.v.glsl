#version 450

in vec4 v_coord;
uniform mat4 m, v, p;

void main()
{
  mat4 mvp = p*v*m;
  gl_Position = mvp * v_coord;
}
