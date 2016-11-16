#ifndef OBJREADER_H
#define OBJREADER_H

#include "SimpleMesh/simplemesh_global.h"
#include "SimpleMesh/mesh.h"
#include <string>

using namespace std;

namespace SimpleMesh{

namespace IO {

class SIMPLEMESH_API OBJReader
{
private:
    string filename;
    bool read_ok;

public:
    OBJReader(const char* file);
    bool is_ok() {return read_ok;}
    void read_mesh(Mesh& mesh);
};

} //namespace IO
} //namespace SimpleMesh
#endif // OBJREADER_H
