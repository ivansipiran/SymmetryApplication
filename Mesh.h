#ifndef __MESH_H_
#define __MESH_H_

/*#include <openvdb/openvdb.h>
#include <openvdb/util/Util.h>
#include <openvdb/tools/MeshToVolume.h>*/
#include "Vertex.h"
#include "Triangle.h"

namespace SymArchData{


class Mesh{
    public:
     typedef enum{
         TRIANGLE_MESH = 0x1,
         POINT_CLOUD = 0x2
     }MeshType;

    private:
        Vertex*                     vertices;
        Triangle*                   triangles;

        float* p_Vertices;
        float* p_Normals;
        unsigned int* p_Elements;

        unsigned int                num_vertices;
        unsigned int                num_triangles;
        Mesh*                       bounding_box;

        float                       rgb[3];
        float                       E[3];

        MeshType                    type_mesh;
        bool                        is_visible;

        //openvdb::FloatGrid::Ptr     grid;
        bool                        hasGrid;

        bool                        plane_visible;
        float err;

        vector<vector<Vertex*> > components;
        vector<unsigned int>         vertex_components;
        vector<bool>                 marked_components;

   public:

        Mesh();
        Mesh(const char* filename, bool hGrid=false);
        ~Mesh();

        float* constVertices();
        float* constNormals();
        unsigned int* constElements();

        void           load_from_file(const char* filename);
        void           load_from_off_file(const char* filename);
        void           load_from_obj_file(const char* filename);

        Vertex*        get_vertices(){return vertices;}
        Triangle*      get_triangles(){return triangles;}

        unsigned int   get_number_vertices(){ return num_vertices;}
        unsigned int   get_number_triangles(){return num_triangles;}
        unsigned int   get_number_fractured_components(){ return components.size();}
        unsigned int   get_vertex_component(int index){return vertex_components[index];}
        void           save_component(char* filename, unsigned int pos);

        void           decompose_connected_components();
        vector<vector<Vertex*> >& get_components(){return components;}
        void           mark_component(unsigned int index){marked_components[index]=true;}
        Mesh*        remove_marked_components();

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

        Vertex*     sample_points(unsigned int num_samples);

        void        save_to_ply(const char* filename);
        void        save_weights(const char* filename);

        void        set_color(float r, float g, float b){rgb[0] = r; rgb[1] = g; rgb[2] = b;}
        float*      get_color(){return rgb;}
        void        set_scale(float ex, float ey, float ez){E[0] = ex; E[1] = ey; E[2] = ez;}
        float*     get_scale(){return E;}

        Mesh*       get_scaled_version();
        void        apply_scale();

        void        set_type(MeshType flag) {type_mesh = flag;}
        MeshType    get_type(){ return type_mesh;}

        void        set_visible(bool t) {is_visible = t;}
        bool        visible(){return is_visible;}

        void        save_to_off(const char* filename, bool flip=false);
        void        save_to_off_with_color(const char* filename);

        float       computeDiagonal();

        //openvdb::FloatGrid::Ptr getGrid(float factor);
        //void        computeGrid(float factor);
        //void        recomputeGrid(float factor);
        void        set_plane_visibility(bool flag){plane_visible=flag;}

        void        set_error(float er){err = er;}
        float      get_error(){return err;}

        float      get_volume();

        void        compute_fractured_components();

        void        flip_faces();
};

}
#endif
