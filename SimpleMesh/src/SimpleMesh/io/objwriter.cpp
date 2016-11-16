#include "SimpleMesh/io/objwriter.h"

namespace SimpleMesh {

namespace IO {

OBJWriter::OBJWriter(const char* file)
{
    filename = file;
}

void OBJWriter :: write_mesh(Mesh &mesh)
{

}

} //namespace IO
} //namespace SimpleMesh

