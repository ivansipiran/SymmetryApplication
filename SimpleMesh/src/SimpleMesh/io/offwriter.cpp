#include <fstream>
#include "SimpleMesh/io/offwriter.h"


namespace SimpleMesh {

namespace IO {

OFFWriter::OFFWriter(const char* file)
{
    filename = string(file);
}

void OFFWriter::write_mesh(Mesh &mesh)
{
    ofstream out(filename.c_str());

    out << "OFF" << endl;
    out << mesh.get_number_vertices() << " " << mesh.get_number_triangles() << " 0" << endl;
    Vertex* vertices = mesh.get_vertices();
    for(unsigned int i = 0; i < mesh.get_number_vertices(); i++){
        out << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << endl;
    }

    Triangle* triangles = mesh.get_triangles();
    for(unsigned int i = 0; i < mesh.get_number_triangles(); i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        out << "3 " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
    }

    out.close();
}

void OFFWriter::write_colored_mesh(Mesh &mesh)
{
    ofstream out(filename.c_str());

    out << "COFF" << endl;
    out << mesh.get_number_vertices() << " " << mesh.get_number_triangles() << " 0" << endl;
    Vertex* vertices = mesh.get_vertices();
    for(unsigned int i = 0; i < mesh.get_number_vertices(); i++){
        unsigned int red = (unsigned int)255*vertices[i].red();
        unsigned int green = (unsigned int)255*vertices[i].green();
        unsigned int blue = (unsigned int)255*vertices[i].blue();
        if(red > 255)   red = 255;
        if(red < 0)     red = 0;
        if(green > 255)   green = 255;
        if(green < 0)     green = 0;
        if(blue > 255)   blue = 255;
        if(blue < 0)     blue = 0;

        out << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << " " << red << " " << green << " " << blue << endl;
    }

    Triangle* triangles = mesh.get_triangles();
    for(unsigned int i = 0; i < mesh.get_number_triangles(); i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        out << "3 " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
    }

    out.close();
}

} // namespace IO
} //namespace SimpleMesh

