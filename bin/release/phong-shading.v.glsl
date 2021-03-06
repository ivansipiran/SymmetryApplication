#version 450

in vec3 v_coord;
in vec3 v_normal;
out vec4 position;  // position of the vertex (and fragment) in world space
out vec3 varyingNormalDirection;  // surface normal vector in world space
uniform mat4 m, v, p;
uniform mat3 m_3x3_inv_transp;


void main()
{
  position = m * vec4(v_coord, 1.0);
  varyingNormalDirection = normalize(m_3x3_inv_transp * v_normal);

  mat4 mvp = p*v*m;
  gl_Position = mvp * vec4(v_coord, 1.0);
}
