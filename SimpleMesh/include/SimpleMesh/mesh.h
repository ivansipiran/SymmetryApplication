#ifndef MESH_H
#define MESH_H

#include "SimpleMesh/simplemesh_global.h"
#include "SimpleMesh/vertex.h"
#include "SimpleMesh/triangle.h"

namespace SimpleMesh{

class SIMPLEMESH_API Mesh
{
protected:
    Vertex*         vertices;
    Triangle*       triangles;
    unsigned int    num_vertices;
    unsigned int    num_triangles;

    float xmin, xmax, ymin, ymax, zmin, zmax;
    float diag;

    bool flag_normals, flag_color;

public:
    Mesh();
    ~Mesh();

    Vertex*     get_vertices() {return vertices;}
    Triangle*   get_triangles() {return triangles;}

    unsigned int get_number_vertices() { return num_vertices; }
    unsigned int get_number_triangles() { return num_triangles;}

    void add_vertex(unsigned int pos, float x, float y, float z);
    void add_triangle(unsigned int pos, unsigned int p1, unsigned int p2, unsigned int p3);
    void add_vertex(unsigned int pos, Vertex v);
    void set_vertex(unsigned int pos, float x, float y, float z);
    void add_normal(unsigned int pos, float nx, float ny, float nz);
    void add_color(unsigned int pos, float r, float g, float b);

    void set_number_vertices(unsigned int nv);
    void set_number_triangles(unsigned int nt);
    void compute_normals();

    float get_diagonal();
    float get_area();

    bool has_normals(){return flag_normals;}
    bool has_color(){return flag_color;}

    void enable_normals(){flag_normals = true;}
    void enable_color(){flag_color = true;}
};

} //namespace SimpleMesh

#endif // MESH_H
