#ifndef MESHT_H
#define MESHT_H

#include <vector>
#include <string>

using namespace std;

class MeshT
{
public:
    MeshT();

    typedef enum {ERROR_OK = 0, ERROR_READ} error_t;

    error_t load_from_file(const char* filename);
    float* constVertices() const {return vertexBuffer.data();}
    float* constNormals() const {return normalBuffer.data();}
    unsigned int* constElements() const {return indexBuffer.data();}
    bool has_color(){return flag_color;}

    int get_number_vertices(){return num_vertices;}
    int get_number_triangles() {return num_triangles;}

    inline float get_x_at(int index) {return vertexBuffer[3 * index];}
    inline float get_y_at(int index) {return vertexBuffer[3 * index + 1];}
    inline float get_z_at(int index) {return vertexBuffer[3 * index + 2];}

    void compute_normals();

private:
    error_t load_from_off_file(const char* filename);
    error_t load_from_obj_file(const char* filename);

private:
    int num_vertices;
    int num_triangles;

    bool flag_color;
    vector<float> vertexBuffer;
    vector<float> normalBuffer;
    vector<float> colorBuffer;
    vector<unsigned int> indexBuffer;

    error_t error;
    string error_message;
};

#endif // MESHT_H
