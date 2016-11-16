#version 430

layout (location = 0) in vec4 v_coord;
layout (location = 1) in vec3 v_normal;
layout (location = 2) in vec3 v_color;

uniform mat4 m, v, p;

out vec3 Color;

void main(void)
{
    Color = v_color;
    mat4 mvp = p*v*m;
    gl_Position = mvp * v_coord;
}
