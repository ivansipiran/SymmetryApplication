/**
 * From the OpenGL Programming wikibook: http://en.wikibooks.org/wiki/OpenGL_Programming
 * This file is in the public domain.
 * Contributors: Martin Kraus, Sylvain Beucler
 */
attribute vec4 v_coord;
uniform mat4 m, v, p;

void main()
{
  mat4 mvp = p*v*m;
  gl_Position = mvp * v_coord;
}
