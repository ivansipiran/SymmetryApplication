#include "mesht.h"

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <sstream>
#include <fstream>

using namespace boost;

MeshT::MeshT() :
    flag_color(false),
    error(0)
{

}

MeshT::error_t MeshT::load_from_file(const char *filename){
    string str(filename);

    if(str.substr(str.length()-3, 3).compare("off")==0 || str.substr(str.length()-3, 3).compare("OFF")==0)
        error = load_from_off_file(filename);
    else if(str.substr(str.length()-3, 3).compare("obj")==0 || str.substr(str.length()-3, 3).compare("OBJ")==0)
        error = load_from_obj_file(filename);
    else{
        error = ERROR_READ;
        error_message = "Unkown file format. Only OFF and OBJ is supported";
    }

    return error;
}

MeshT::error_t MeshT :: load_from_obj_file (const char* filename){
        ifstream in(filename);

        string str;

        //vector<Vertex> vert;
        //vector<Triangle> fac;

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

            if (tok[0].compare("v") == 0){
                //Vertex v;
                //v.setX(atof(tok[1].c_str()));
                //v.setY(atof(tok[2].c_str()));
                //v.setZ(atof(tok[3].c_str()));
                //vert.push_back(v);
                //cout << v.x() << " " << v.y() << " " << v.z() << endl;
                vertexBuffer.push_back(atof(tok[1].c_str()));
                vertexBuffer.push_back(atof(tok[2].c_str()));
                vertexBuffer.push_back(atof(tok[3].c_str()));
            }
            else if (tok[0].compare("f") == 0){
                int numVertex = tok.size() - 1;
                if (numVertex < 3){
                    error = ERROR_READ;
                    stringstream ss;
                    ss << "Error while reading face on line:" << lineNumber;
                    error_message = ss.str();
                    return error;
                }
                vector<unsigned int> indices;
                char_separator<char> sep2("/");
                for (int i = 1; i < tok.size(); i++){
                    tokenizer< char_separator<char> > tokens2(tok[i], sep2);
                    tokenizer< char_separator<char> >::iterator beg = tokens2.begin();
                    unsigned int index = (unsigned int)atoi(beg->c_str());

                    vector<unsigned int>::iterator it = find(indices.begin(), indices.end(), index);
                    if (it == indices.end())
                        //indices.push_back(index);
                        indexBuffer.push_back(index-1);
                    else{
                        cout << "Repeating vertices in face list on line:" << lineNumber << "(" << str << ")" << endl;
                    }
                    //cout << index << " ";
                }

                //cout << endl;
                //Triangle f;
                //f.add_vertex(indices[0] - 1);
                //f.add_vertex(indices[1] - 1);
                //f.add_vertex(indices[2] - 1);
                //fac.push_back(f);
            }
            lineNumber++;
        }

        num_vertices = vertexBuffer.size()/3;
        num_triangles = indexBuffer.size()/3;

        //vertices = new Vertex[num_vertices];
        //triangles = new Triangle[num_triangles];

        //for (unsigned int i = 0; i < vert.size(); i++){
        //    vertices[i].setX(vert[i].x());
        //    vertices[i].setY(vert[i].y());
        //    vertices[i].setZ(vert[i].z());

        //    vertices[i].set_index(i);
        //}

        //for (unsigned int i = 0; i < fac.size(); i++){
        //    unsigned int p1, p2, p3;
        //    p1 = fac[i].get_vertex_at(0);
        //    p2 = fac[i].get_vertex_at(1);
        //    p3 = fac[i].get_vertex_at(2);

        //    triangles[i].add_vertex(p1);		triangles[i].add_vertex(p2);		triangles[i].add_vertex(p3);

        //    vertices[p1].add_face(i);	vertices[p2].add_face(i);	vertices[p3].add_face(i);
        //    vertices[p1].addVertex(p2);	vertices[p1].addVertex(p3);
        //    vertices[p2].addVertex(p1);	vertices[p2].addVertex(p3);
        //    vertices[p3].addVertex(p1);	vertices[p3].addVertex(p2);
        //}

        in.close();
        return ERROR_OK;
}

MeshT::error_t MeshT :: load_from_off_file (const char* filename){
    bool flagColor=false, flagNormal=false, flagTexture=false;

    ifstream in(filename);
    string str;

    //Read header
    getline(in, str);
    if(str.find("OFF") == string::npos){
        error = ERROR_READ;
        error_message = "OFF format requires the mark OFF as first line in the file";
        return error;
    }


    if(str.find("ST")!=string::npos)
        flagTexture = true;

    if(str.find("C")!=string::npos)
        flagColor = true;

    if(str.find("N")!=string::npos)
        flagNormal = true;

    //Read header - info mesh
    unsigned int nVertices, nTriangles, nEdges;
    getline(in, str);

    vector<string> tokIni;
    char_separator<char> sepIni(" ");
    tokenizer< char_separator<char> > tokensIni(str, sepIni);
    BOOST_FOREACH(const string& t, tokensIni){
        tokIni.push_back(t);
    }

    assert(tokIni.size() == 3);
    num_vertices = (unsigned int)atoi(tokIni[0].c_str());
    num_triangles = (unsigned int)atoi(tokIni[1].c_str());
    nEdges = (unsigned int)atoi(tokIni[2].c_str());

    //Allocate storage for the mesh
    //set_num_vertices(nVertices);
    //set_num_triangles(nTriangles);

    //Actually reading the mesh structure
    for(unsigned int i = 0; i < num_vertices; i++){
        getline(in, str);
        vector<string> tok;
        char_separator<char> sep(" ");
        tokenizer< char_separator<char> > tokens(str, sep);
        BOOST_FOREACH(const string& t, tokens){
            tok.push_back(t);
        }

        //add_vertex(i, atof(tok[0].c_str()), atof(tok[1].c_str()), atof(tok[2].c_str()));
        vertexBuffer.push_back(atof(tok[0].c_str()));
        vertexBuffer.push_back(atof(tok[1].c_str()));
        vertexBuffer.push_back(atof(tok[2].c_str()));

        if(flagNormal){
            //add_normal(i, atof(tok[3].c_str()), atof(tok[4].c_str()), atof(tok[5].c_str()));
            normalBuffer.push_back(atof(tok[3].c_str()));
            normalBuffer.push_back(atof(tok[4].c_str()));
            normalBuffer.push_back(atof(tok[5].c_str()));
        }
        if(flagColor){
            //add_color(i, (unsigned int)atoi(tok[6].c_str()), (unsigned int)atoi(tok[7].c_str()), (unsigned int)atoi(tok[8].c_str()));
            colorBuffer.push_back(atof(tok[6].c_str())/255);
            colorBuffer.push_back(atof(tok[7].c_str())/255);
            colorBuffer.push_back(atof(tok[8].c_str())/255);
        }
    }
    flag_color = flagColor;

    for(unsigned int i = 0; i < num_triangles; i++){
        getline(in, str);
        vector<string> tok;
        char_separator<char> sep(" ");
        tokenizer< char_separator<char> > tokens(str, sep);
        BOOST_FOREACH(const string& t, tokens){
            tok.push_back(t);
        }

        if(atoi(tok[0].c_str()) != 3){
            error = ERROR_READ;
            stringstream ss;
            ss << "There is no a triangle in line " << i;
            error_message = ss.str();
            return error;
        }

        unsigned int p1, p2, p3;
        p1 = atoi(tok[1].c_str());
        p2 = atoi(tok[2].c_str());
        p3 = atoi(tok[3].c_str());
        //add_triangle(i, p1, p2, p3);
        indexBuffer.push_back(p1);
        indexBuffer.push_back(p2);
        indexBuffer.push_back(p3);
    }

    in.close();
    return ERROR_OK;
}

void MeshT :: compute_normals(){
    vector<float> normal_t;
    vector< vector<unsigned int> > adjacency;

    adjacency.reserve(num_vertices);

    //Compute normals for faces
    for(unsigned int i = 0; i < num_triangles; i++){
        float nx = 0.0, ny = 0.0, nz = 0.0;
        vector<unsigned int> vert;
        vert.push_back(indexBuffer[3 * i]);
        vert.push_back(indexBuffer[3 * i + 1]);
        vert.push_back(indexBuffer[3 * i + 2]);

        adjacency[indexBuffer[3 * i]].push_back(i);
        adjacency[indexBuffer[3 * i] + 1].push_back(i);
        adjacency[indexBuffer[3 * i] + 2].push_back(i);

        unsigned int ind1 = 2;
        for(unsigned int j = 0; j < 3; j++){
            int ind2 = j;
            nx += (get_y_at(vert[ind1]) - get_y_at(vert[ind2]))*(get_z_at(vert[ind1]) + get_z_at(vert[ind2]));
            ny += (get_z_at(vert[ind1]) - get_z_at(vert[ind2]))*(get_x_at(vert[ind1]) + get_x_at(vert[ind2]));
            nz += (get_x_at(vert[ind1]) - get_x_at(vert[ind2]))*(get_y_at(vert[ind1]) + get_y_at(vert[ind2]));
            ind1 = ind2;
        }

        float magnitude = sqrt(nx*nx + ny*ny + nz*nz);
        if(magnitude > 1.0E-6){
            nx /= magnitude;
            ny /= magnitude;
            nz /= magnitude;
        }
        normal_t.push_back(nx);
        normal_t.push_back(ny);
        normal_t.push_back(nz);
    }

    //Compute normals for vertices
    for(unsigned int i = 0; i < num_vertices; i++){
        //vector<unsigned int> fac = adjacency[i];
        float nx = 0.0, ny = 0.0, nz = 0.0;
        for(unsigned int j = 0; j < adjacency[i].size(); j++){
            nx += normal_t[3 * adjacency[i][j]];
            ny += normal_t[3 * adjacency[i][j] + 1];
            nz += normal_t[3 * adjacency[i][j] + 2];;
        }

        nx /= adjacency[i].size();
        ny /= adjacency[i].size();
        nz /= adjacency[i].size();

        float magnitude = sqrt(nx*nx + ny*ny + nz*nz);
        if(magnitude > 1.0E-6){
            nx /= magnitude;
            ny /= magnitude;
            nz /= magnitude;
        }

        //cout << nx << " " << ny <<  " " << nz << endl;
        normalBuffer.push_back(nx);
        normalBuffer.push_back(ny);
        normalBuffer.push_back(nz);
    }
}
