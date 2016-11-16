#include <fstream>
#include <iostream>
#include <tuple>

#include "SimpleMesh/io/objreader.h"
#include "SimpleMesh/vertex.h"
#include "SimpleMesh/triangle.h"
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace boost;
using namespace std;

namespace SimpleMesh{

namespace IO {

OBJReader::OBJReader(const char* file)
{
    read_ok = false;
    ifstream in(file);

    if(!in.is_open()){
        return;
    }

    read_ok = true;
    in.close();
    filename = file;
}

void OBJReader::read_mesh(Mesh& mesh)
{
    ifstream in(filename.c_str());
    string str;

    vector<Vertex> vert;
    vector<Triangle> fac;
    vector< std::tuple<float, float, float> > normals;

    bool flagNormals=false, flagColor=false;

    int lineNumber = 1;
    while (!in.eof()){
         getline(in, str);

         if (str.size() == 0 || str[0] == '#')
            continue;

         //Extraer los tokens y guardarlos en un vector
         vector<string> tok;
         char_separator<char> sep(" ");
         tokenizer< char_separator<char> > tokens(str, sep);
         BOOST_FOREACH(const string& t, tokens){
             tok.push_back(t);
         }

         if (tok[tok.size() - 1].find('\\') != string::npos){
             string str1;
             lineNumber++;
             getline(in, str1);
             tokenizer< char_separator<char> > tokens3(str1, sep);
             BOOST_FOREACH(const string& t, tokens3){
                 tok.push_back(t);
             }
         }

         if (tok[0].compare("vn") == 0){
             float nx = atof(tok[1].c_str());
             float ny = atof(tok[2].c_str());
             float nz = atof(tok[3].c_str());

             float mag = sqrt(nx*nx + ny*ny + nz*nz);
             normals.push_back(std::make_tuple(nx/mag, ny/mag, nz/mag));
             flagNormals = true;
         }else if (tok[0].compare("v") == 0){
             Vertex v;
             v.setX(atof(tok[1].c_str()));
             v.setY(atof(tok[2].c_str()));
             v.setZ(atof(tok[3].c_str()));
             if(tok.size()==7){ //Vertex contains color
                 v.set_color(atof(tok[4].c_str()), atof(tok[5].c_str()), atof(tok[6].c_str()));
                 flagColor=true;
             }
             vert.push_back(v);
         }
         else if (tok[0].compare("f") == 0){
             int numVertex = tok.size() - 1;
             if (numVertex < 3){
                 cout << "Error while reading face on line:" << lineNumber << endl;
                 exit(EXIT_FAILURE);
             }
             vector<unsigned int> indices;
             char_separator<char> sep2("/");
             for (int i = 1; i < tok.size(); i++){
                tokenizer< char_separator<char> > tokens2(tok[i], sep2);
                tokenizer< char_separator<char> >::iterator beg = tokens2.begin();
                unsigned int index = (unsigned int)atoi(beg->c_str());

                vector<unsigned int>::iterator it = find(indices.begin(), indices.end(), index);
                if (it == indices.end())
                    indices.push_back(index);
                else{
                    cout << "Repeating vertices in face list on line:" << lineNumber << "(" << str << ")" << endl;
                }
              }

              Triangle f;
              f.add_vertex(indices[0] - 1);
              f.add_vertex(indices[1] - 1);
              f.add_vertex(indices[2] - 1);
              fac.push_back(f);
          }
          lineNumber++;
      }

      mesh.set_number_vertices(vert.size());
      mesh.set_number_triangles(fac.size());

      for (unsigned int i = 0; i < vert.size(); i++){
          mesh.add_vertex(i, vert[i].x(), vert[i].y(), vert[i].z());
          if(flagNormals)
            mesh.add_normal(i, std::get<0>(normals[i]), std::get<1>(normals[i]), std::get<2>(normals[i]));
          if(flagColor)
            mesh.add_color(i, vert[i].red(), vert[i].green(), vert[i].blue());

      }

      for (unsigned int i = 0; i < fac.size(); i++){
          mesh.add_triangle(i, fac[i].get_vertex_at(0), fac[i].get_vertex_at(1), fac[i].get_vertex_at(2));
      }

      if(flagNormals)
          mesh.enable_normals();
      if(flagColor)
          mesh.enable_color();

      in.close();
}

} //namespace IO
}//namespace SimpleMesh
