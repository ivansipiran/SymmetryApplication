#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "SimpleMesh/simplemesh_global.h"
#include <vector>

using namespace std;

namespace SimpleMesh{

class SIMPLEMESH_API Triangle
{
protected:
    vector<unsigned int >   vertices;
    float normal[3];

public:
    Triangle(){ normal[0] = normal[1] = normal[2] = 0.0;}

    void add_vertex(unsigned int v) {vertices.push_back(v);}
    vector<unsigned int>& get_vertices() {return vertices;}
    unsigned int get_vertex_at(unsigned int pos) {return vertices[pos];}

    void set_normal(float nx, float ny, float nz) {normal[0] = nx; normal[1] = ny; normal[2] = nz;}

    float* get_normal() {return normal;}
    float nx() const {return normal[0];}
    float ny() const {return normal[1];}
    float nz() const {return normal[2];}

    float& nx() {return normal[0];}
    float& ny() {return normal[1];}
    float& nz() {return normal[2];}
};
} //namespace SimpleMesh
#endif // TRIANGLE_H
