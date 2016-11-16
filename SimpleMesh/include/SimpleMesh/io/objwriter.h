#ifndef OBJWRITER_H
#define OBJWRITER_H

#include "SimpleMesh/simplemesh_global.h"
#include <string>
#include "SimpleMesh/mesh.h"

using namespace std;

namespace SimpleMesh {

namespace IO {

class SIMPLEMESH_API OBJWriter
{
    string filename;
public:
    OBJWriter(const char* file);
    void write_mesh(Mesh& mesh);
};

} //namespace IO
} //namespace SimpleMesh
#endif // OBJWRITER_H
