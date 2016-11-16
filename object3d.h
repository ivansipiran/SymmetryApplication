#ifndef OBJECT3D_H
#define OBJECT3D_H

#include <cstddef>
#include "Vertex.h"
#include "Triangle.h"
#include "Mesh.h"

/**
 * The class Object3D stores the basic information for a 3D object. Basic functionalities include:
 *
 *   - Creation sub-routines for vertices, triangles, normals and colors
 *   - Compute normals given the vertices and the conectivity
 */

namespace SymArchData{

class Object3D
{
    public:
        typedef enum{ TRIANGLE_MESH = 0x1, POINT_CLOUD = 0x2} MeshType;

    private:
        Vertex*                     vertices;
        Triangle*                   triangles;

        unsigned int                num_vertices;
        unsigned int                num_triangles;
        Mesh*                       bounding_box;

        float                       rgb[3];
        float                       E[3];

        MeshType                    type_mesh;
        bool                        is_visible;

        float err;

    public:
        Object3D();
        ~Object3D();

        Vertex*        get_vertices(){return vertices;}
        Triangle*      get_triangles(){return triangles;}

        unsigned int   get_number_vertices(){ return num_vertices;}
        unsigned int   get_number_triangles(){return num_triangles;}

        void        add_vertex(unsigned int pos, float x, float y, float z);
        void        add_triangle(unsigned int pos, unsigned int p1, unsigned int p2, unsigned int p3);
        void        add_vertex(Vertex v, unsigned int pos);
        void        set_vertex(unsigned int pos, float x, float y, float z);
        void        add_normal(unsigned int pos, float nx, float ny, float nz);
        void        add_color(unsigned int pos, unsigned int r, unsigned int g, unsigned int b);
        void        add_interest(unsigned int pos, bool in){ vertices[pos].setInterest(in);}

        void        set_num_vertices(unsigned int nv);
        void        set_num_triangles(unsigned int nt);

        void        compute_normals();
        float       normalize_area();

        void        set_color(float r, float g, float b){rgb[0] = r; rgb[1] = g; rgb[2] = b;}
        float*      get_color(){return rgb;}
        void        set_scale(float ex, float ey, float ez){E[0] = ex; E[1] = ey; E[2] = ez;}
        float*      get_scale(){return E;}

        void        set_type(MeshType flag) {type_mesh = flag;}
        MeshType    get_type(){ return type_mesh;}

        float       computeDiagonal();
};

}
#endif // OBJECT3D_H
