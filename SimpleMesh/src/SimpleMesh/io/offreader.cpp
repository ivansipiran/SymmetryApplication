#include <iostream>
#include <fstream>

#include "SimpleMesh/io/offreader.h"
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

namespace SimpleMesh{

namespace IO {

OFFReader::OFFReader(const char* file)
{
    read_ok = false;
    ifstream in(file);

    if(!in.is_open()){
        return;
    }
	
	cout << "File " << file << " - opened!" << endl;

    string str;

    getline(in, str);
    if(str.find("OFF") != string::npos){
        read_ok = true;
        filename = file;
	}
	else{
		cout << "File " << file << " format! - Bad" << endl;
	}

	cout << "File " << file << " - " << read_ok << endl;
    in.close();
}

void OFFReader::read_mesh(Mesh& mesh)
{
    bool flagColor=false, flagNormal=false, flagTexture=false;

    ifstream in(filename.c_str());
    string str;

    //Read header
    getline(in, str);
    #ifdef DEBUG
        cout << str << endl;
    #endif // DEBUG

    if(str.find("ST")!=string::npos)
        flagTexture = true;

    if(str.find("C")!=string::npos)
        flagColor = true;

    if(str.find("N")!=string::npos)
        flagNormal = true;

    //Read header - info mesh
    unsigned int nVertices, nTriangles, nEdges;
    getline(in, str);
    #ifdef DEBUG
        cout << str << endl;
    #endif // DEBUG

    vector<string> tokIni;
    char_separator<char> sepIni(" ");
    tokenizer< char_separator<char> > tokensIni(str, sepIni);
    BOOST_FOREACH(const string& t, tokensIni){
        tokIni.push_back(t);
    }

    assert(tokIni.size() == 3);
    nVertices = (unsigned int)atoi(tokIni[0].c_str());
    nTriangles = (unsigned int)atoi(tokIni[1].c_str());
    nEdges = (unsigned int)atoi(tokIni[2].c_str());

    #ifdef DEBUG
        cout << nVertices << " - " << nTriangles << " - " << nEdges << endl;
    #endif // DEBUG

    //Allocate storage for the mesh
    mesh.set_number_vertices(nVertices);
    mesh.set_number_triangles(nTriangles);

    #ifdef DEBUG
	cout << "#vertex:" << mesh.get_number_vertices() << " - set!" << endl;
	cout << "#triangle:" << mesh.get_number_triangles() << " - set!" << endl;
    #endif // DEBUG

    //Actually reading the mesh structure
    for(unsigned int i = 0; i < mesh.get_number_vertices(); i++){
        getline(in, str);
        vector<string> tok;
        char_separator<char> sep(" ");
        tokenizer< char_separator<char> > tokens(str, sep);
        BOOST_FOREACH(const string& t, tokens){
            tok.push_back(t);
        }

        mesh.add_vertex(i, atof(tok[0].c_str()), atof(tok[1].c_str()), atof(tok[2].c_str()));

        if(flagNormal)
            mesh.add_normal(i, atof(tok[3].c_str()), atof(tok[4].c_str()), atof(tok[5].c_str()));
        if(flagColor){
            if(!flagNormal){
                mesh.add_color(i, atof(tok[3].c_str())/255.0, atof(tok[4].c_str())/255.0, atof(tok[5].c_str())/255.0);
            }
            if(flagNormal){
                mesh.add_color(i, atof(tok[6].c_str())/255.0, atof(tok[7].c_str())/255.0, atof(tok[8].c_str())/255.0);
            }
        }
    }

    for(unsigned int i = 0; i < mesh.get_number_triangles(); i++){
        getline(in, str);
        vector<string> tok;
        char_separator<char> sep(" ");
        tokenizer< char_separator<char> > tokens(str, sep);
        BOOST_FOREACH(const string& t, tokens){
            tok.push_back(t);
        }

        assert(atoi(tok[0].c_str()) == 3);
        unsigned int p1, p2, p3;
        p1 = atoi(tok[1].c_str());
        p2 = atoi(tok[2].c_str());
        p3 = atoi(tok[3].c_str());
        mesh.add_triangle(i, p1, p2, p3);

    }

    if(flagNormal)
        mesh.enable_normals();
    if(flagColor)
        mesh.enable_color();

    in.close();
}

} //namespace IO
}//namespace SimpleMesh
