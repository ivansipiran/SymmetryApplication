#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <map>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include "Mesh.h"
#include "util.h"

using namespace std;
using namespace boost;

namespace SymArchData{


typedef Vertex Vector;

float determinant(Vector a, Vector b){
	return a.x()*(b.y() - b.z()) - a.y()*(b.x() - b.z()) + a.z()*(b.x() - b.y());
}

Mesh :: Mesh (){
    vertices = NULL;
    triangles = NULL;
    num_vertices = 0;
    num_triangles = 0;
    bounding_box = NULL;
    type_mesh = TRIANGLE_MESH;
    is_visible = true;
    E[0] = E[1] = E[2] = 1.0f;
    hasGrid = false;
    plane_visible = false;
}

Mesh :: Mesh(const char* filename, bool hGrid){
    vertices = NULL;
    triangles = NULL;
	num_vertices = 0;
	num_triangles = 0;
	bounding_box = NULL;
    type_mesh = TRIANGLE_MESH;

	load_from_file(filename);

    compute_normals();
    E[0] = E[1] = E[2] = 1.0;
    is_visible = true;
    hasGrid = hGrid;
    plane_visible = false;

    if(hasGrid){

    }
}

/*void Mesh :: computeGrid(float factor){

    openvdb::math::Transform::Ptr t = openvdb::math::Transform::createLinearTransform(factor);

    std::vector<openvdb::Vec3s> pointList;
    std::vector<openvdb::Vec4I> polygonList;

    for(int i = 0; i < get_number_vertices(); i++){
        openvdb::Vec3s xyz(0.0, 0.0, 0.0);
        xyz = t->worldToIndex(openvdb::Vec3s(vertices[i].x(), vertices[i].y(), vertices[i].z()));
        pointList.push_back(xyz);
    }

    for(int i = 0; i < get_number_triangles(); i++){
        std::vector<unsigned int> vert = triangles[i].get_vertices();
        polygonList.push_back(openvdb::Vec4I(vert[0], vert[1], vert[2], openvdb::util::INVALID_IDX));
    }

    grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*t, pointList, polygonList);
}

openvdb::FloatGrid::Ptr Mesh :: getGrid(float factor){
    if(!hasGrid){
        computeGrid(factor);
    }
    return grid;
}

void Mesh :: recomputeGrid(float factor){

    computeGrid(factor);
}*/

Mesh :: ~Mesh(){
    if(triangles){
        delete[] triangles;
        triangles = NULL;
    }

    if(vertices){
        delete[] vertices;
        vertices = NULL;
    }
    if(bounding_box)
        delete bounding_box;

    if(p_Vertices)
        delete p_Vertices;

    if(p_Normals)
        delete p_Normals;

    if(p_Elements)
        delete p_Elements;
}

void Mesh :: load_from_file(const char* filename){
    string str(filename);

    if(str.substr(str.length()-3, 3).compare("off")==0 || str.substr(str.length()-3, 3).compare("OFF")==0)
        load_from_off_file(filename);
    else if(str.substr(str.length()-3, 3).compare("obj")==0 || str.substr(str.length()-3, 3).compare("OBJ")==0)
        load_from_obj_file(filename);
    else{
        cout << "Unknown format" << endl;
        exit(EXIT_FAILURE);
    }

}

void Mesh :: load_from_obj_file (const char* filename){
        ifstream in(filename);

        string str;

        vector<Vertex> vert;
        vector<Triangle> fac;

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
                Vertex v;
                v.setX(atof(tok[1].c_str()));
                v.setY(atof(tok[2].c_str()));
                v.setZ(atof(tok[3].c_str()));
                vert.push_back(v);
                //cout << v.x() << " " << v.y() << " " << v.z() << endl;
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
                    //cout << index << " ";
                }

                //cout << endl;
                Triangle f;
                f.add_vertex(indices[0] - 1);
                f.add_vertex(indices[1] - 1);
                f.add_vertex(indices[2] - 1);
                fac.push_back(f);
            }
            lineNumber++;
        }

        num_vertices = vert.size();
        num_triangles = fac.size();

        vertices = new Vertex[num_vertices];
        triangles = new Triangle[num_triangles];

        for (unsigned int i = 0; i < vert.size(); i++){
            vertices[i].setX(vert[i].x());
            vertices[i].setY(vert[i].y());
            vertices[i].setZ(vert[i].z());

            vertices[i].set_index(i);
        }

        for (unsigned int i = 0; i < fac.size(); i++){
            unsigned int p1, p2, p3;
            p1 = fac[i].get_vertex_at(0);
            p2 = fac[i].get_vertex_at(1);
            p3 = fac[i].get_vertex_at(2);

            triangles[i].add_vertex(p1);		triangles[i].add_vertex(p2);		triangles[i].add_vertex(p3);

            vertices[p1].add_face(i);	vertices[p2].add_face(i);	vertices[p3].add_face(i);
            vertices[p1].addVertex(p2);	vertices[p1].addVertex(p3);
            vertices[p2].addVertex(p1);	vertices[p2].addVertex(p3);
            vertices[p3].addVertex(p1);	vertices[p3].addVertex(p2);
        }

        in.close();
}

void Mesh :: load_from_off_file (const char* filename){
    bool flagColor=false, flagNormal=false, flagTexture=false;

    ifstream in(filename);
    string str;

    //Read header
    getline(in, str);
    assert(str.find("OFF") != string::npos);

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
    nVertices = (unsigned int)atoi(tokIni[0].c_str());
    nTriangles = (unsigned int)atoi(tokIni[1].c_str());
    nEdges = (unsigned int)atoi(tokIni[2].c_str());

    //Allocate storage for the mesh
    set_num_vertices(nVertices);
    set_num_triangles(nTriangles);

    //Actually reading the mesh structure
    for(unsigned int i = 0; i < num_vertices; i++){
        getline(in, str);
        vector<string> tok;
        char_separator<char> sep(" ");
        tokenizer< char_separator<char> > tokens(str, sep);
        BOOST_FOREACH(const string& t, tokens){
            tok.push_back(t);
        }

        add_vertex(i, atof(tok[0].c_str()), atof(tok[1].c_str()), atof(tok[2].c_str()));

        if(flagNormal)
            add_normal(i, atof(tok[3].c_str()), atof(tok[4].c_str()), atof(tok[5].c_str()));
        if(flagColor)
            add_color(i, (unsigned int)atoi(tok[6].c_str()), (unsigned int)atoi(tok[7].c_str()), (unsigned int)atoi(tok[8].c_str()));
    }

    for(unsigned int i = 0; i < num_triangles; i++){
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
        add_triangle(i, p1, p2, p3);
    }

    in.close();
}

void Mesh :: compute_normals(){
    //Compute normals for faces
    for(unsigned int i = 0; i < num_triangles; i++){
        float nx = 0.0, ny = 0.0, nz = 0.0;
        vector<unsigned int> vert = triangles[i].get_vertices();
        unsigned int ind1 = 2;
        for(unsigned int j = 0; j < 3; j++){
            int ind2 = j;
            nx += (vertices[vert[ind1]].y() - vertices[vert[ind2]].y())*(vertices[vert[ind1]].z() + vertices[vert[ind2]].z());
            ny += (vertices[vert[ind1]].z() - vertices[vert[ind2]].z())*(vertices[vert[ind1]].x() + vertices[vert[ind2]].x());
            nz += (vertices[vert[ind1]].x() - vertices[vert[ind2]].x())*(vertices[vert[ind1]].y() + vertices[vert[ind2]].y());
            ind1 = ind2;
        }

        float magnitude = sqrt(nx*nx + ny*ny + nz*nz);
        if(magnitude > 1.0E-6){
            nx /= magnitude;
            ny /= magnitude;
            nz /= magnitude;
        }
        triangles[i].set_normal(nx, ny, nz);
    }

    //Compute normals for vertices
    for(unsigned int i = 0; i < num_vertices; i++){
        vector<unsigned int> fac = vertices[i].get_faces();
        float nx = 0.0, ny = 0.0, nz = 0.0;
        for(unsigned int j = 0; j < fac.size(); j++){
            nx += triangles[fac[j]].nx();
            ny += triangles[fac[j]].ny();
            nz += triangles[fac[j]].nz();
        }

        nx /= fac.size();
        ny /= fac.size();
        nz /= fac.size();

        float magnitude = sqrt(nx*nx + ny*ny + nz*nz);
        if(magnitude > 1.0E-6){
            nx /= magnitude;
            ny /= magnitude;
            nz /= magnitude;
        }

        vertices[i].set_normal(nx, ny, nz);
    }
}

void Mesh :: set_num_vertices(unsigned int nv){
    num_vertices = nv;
    vertices = new Vertex[num_vertices];
}

void Mesh :: set_num_triangles(unsigned int nt){
    num_triangles = nt;
    triangles = new Triangle[num_triangles];
}

void Mesh ::add_vertex(unsigned int pos, float x, float y, float z){
    vertices[pos].x() = x;
    vertices[pos].y() = y;
    vertices[pos].z() = z;
    vertices[pos].set_index(pos);
}

void Mesh :: set_vertex(unsigned int pos, float x, float y, float z){
    vertices[pos].x() = x;
    vertices[pos].y() = y;
    vertices[pos].z() = z;
}

void Mesh :: add_vertex(Vertex v, unsigned int pos)
{
    vertices[pos].x() = v.x();
    vertices[pos].y() = v.y();
    vertices[pos].z() = v.z();
    vertices[pos].set_index(pos);
    vertices[pos].set_weight(v.get_weight());
    vertices[pos].compute_color();
}

void Mesh :: add_triangle(unsigned int pos, unsigned int p1, unsigned int p2, unsigned int p3){
    triangles[pos].add_vertex(p1);
    triangles[pos].add_vertex(p2);
    triangles[pos].add_vertex(p3);

    vertices[p1].add_face(pos);
    vertices[p2].add_face(pos);
    vertices[p3].add_face(pos);

    vertices[p1].addVertex(p2);	vertices[p1].addVertex(p3);
    vertices[p2].addVertex(p1);	vertices[p2].addVertex(p3);
    vertices[p3].addVertex(p1);	vertices[p3].addVertex(p2);
}

void Mesh :: add_normal(unsigned int pos, float nx, float ny, float nz){
    vertices[pos].set_normal(nx, ny, nz);
}

void Mesh :: add_color(unsigned int pos, unsigned int r, unsigned int g, unsigned int b){
    vertices[pos].set_color(r, g, b);
}

float Mesh::normalize_area(){
    float sum = 0.0;
    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        Vector XV(vertices[vert[0]].x(), vertices[vert[1]].x(), vertices[vert[2]].x());
        Vector YV(vertices[vert[0]].y(), vertices[vert[1]].y(), vertices[vert[2]].y());
        Vector ZV(vertices[vert[0]].z(), vertices[vert[1]].z(), vertices[vert[2]].z());

        float det1 = determinant(XV, YV);
        float det2 = determinant(YV, ZV);
        float det3 = determinant(ZV, XV);

        float area = 0.5*sqrt(det1*det1 + det2*det2 + det3*det3);
        sum = sum + area;
    }

    float sc = 1.0/sqrt(sum);
    E[0] = E[1] = E[2] = sc;
    apply_scale();
    return sc;
}

float Mesh :: get_volume(){
    float sum = 0.0;

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        float v321 = vertices[vert[2]].x() * vertices[vert[1]].y() * vertices[vert[0]].z();
        float v231 = vertices[vert[1]].x() * vertices[vert[2]].y() * vertices[vert[0]].z();
        float v312 = vertices[vert[2]].x() * vertices[vert[0]].y() * vertices[vert[1]].z();
        float v132 = vertices[vert[0]].x() * vertices[vert[2]].y() * vertices[vert[1]].z();
        float v213 = vertices[vert[1]].x() * vertices[vert[0]].y() * vertices[vert[2]].z();
        float v123 = vertices[vert[0]].x() * vertices[vert[1]].y() * vertices[vert[2]].z();
        sum += (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
    }

    return sum;
}

Vertex* Mesh :: sample_points(unsigned int num_samples){
    compute_normals();

    //std::cout << "Normals computed" << std::endl;

    float* accum = new float[num_triangles];
    float sum = 0.0;

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        Vector XV(vertices[vert[0]].x(), vertices[vert[1]].x(), vertices[vert[2]].x());
		Vector YV(vertices[vert[0]].y(), vertices[vert[1]].y(), vertices[vert[2]].y());
		Vector ZV(vertices[vert[0]].z(), vertices[vert[1]].z(), vertices[vert[2]].z());

        float det1 = determinant(XV, YV);
        float det2 = determinant(YV, ZV);
        float det3 = determinant(ZV, XV);

        float area = 0.5f*sqrt(det1*det1 + det2*det2 + det3*det3);
		accum[i] = sum + area;
		sum = accum[i];
    }

    //std::cout << "Accumulator computed" << std::endl;

    unsigned int curr_num_samples = 0;
    float max_accum = accum[num_triangles - 1];

    Vertex* point_cloud = new Vertex[num_samples];

    srand(time(NULL));

    //std::cout << "Start sampling" << std::endl;
    while(curr_num_samples < num_samples){
        float value = (rand()*max_accum)/RAND_MAX;
        unsigned int tt = Util::binary_search(accum, num_triangles, value);

        vector<unsigned int> vert1 = triangles[tt].get_vertices();

        float r1 = ((float)rand())/RAND_MAX;
        float r2 = ((float)rand())/RAND_MAX;

        float sqr = sqrt(r1);
        float factor1 = (1 - sqr);
        float factor2 = sqr * (1 - r2);
        float factor3 = sqr * r2;

        float x1 =  factor1 * vertices[vert1[0]].x() + factor2 * vertices[vert1[1]].x() + factor3 * vertices[vert1[2]].x();
        float y1 =  factor1 * vertices[vert1[0]].y() + factor2 * vertices[vert1[1]].y() + factor3 * vertices[vert1[2]].y();
        float z1 =  factor1 * vertices[vert1[0]].z() + factor2 * vertices[vert1[1]].z() + factor3 * vertices[vert1[2]].z();

		point_cloud[curr_num_samples].x() = x1;
		point_cloud[curr_num_samples].y() = y1;
		point_cloud[curr_num_samples].z() = z1;

		point_cloud[curr_num_samples].nx() = triangles[tt].nx();
		point_cloud[curr_num_samples].ny() = triangles[tt].ny();
		point_cloud[curr_num_samples].nz() = triangles[tt].nz();

		curr_num_samples++;
    }
    //std::cout << "Finish sampling" << std::endl;

    delete[] accum;
    return point_cloud;
}

void Mesh :: save_to_ply(const char* filename){
    ofstream in(filename);

    in << "ply" << endl;
    in << "format ascii 1.0" << endl;
    in << "element vertex " << num_vertices << endl;
    in << "property float x" << endl;
    in << "property float y" << endl;
    in << "property float z" << endl;
    in << "property uint red" << endl;
    in << "property uint green" << endl;
    in << "property uint blue" << endl;
    in << "element face " << num_triangles << endl;
    in << "property list uchar int vertex_indices" << endl;
    in << "end_header" << endl;

    for(unsigned int i = 0; i < num_vertices; i++){
        in << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << " " << vertices[i].red() << " " << vertices[i].green() << " " << vertices[i].blue() << endl;
    }

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        in << "3 " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
    }

    in.close();
}

void Mesh :: save_weights(const char* filename){
    ofstream out(filename);
    for(int i = 0; i < num_vertices; i++)
        out << vertices[i].get_weight() << endl;
    out.close();
}

Mesh* Mesh::get_scaled_version()
{
    Mesh* new_mesh = new Mesh();

    new_mesh->set_num_vertices(num_vertices);
    new_mesh->set_num_triangles(num_triangles);

    for(unsigned int i = 0; i < num_vertices; i++){
        new_mesh->add_vertex(E[0]*vertices[i].x(), E[1]*vertices[i].y(), E[2]*vertices[i].z(), i);
    }

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        new_mesh->add_triangle(i, vert[0], vert[1], vert[2]);
    }
    new_mesh->set_type(type_mesh);

    return new_mesh;
}

void Mesh :: flip_faces(){
    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        unsigned int a = vert[1];
        unsigned int b = vert[2];
        vert[1] = b;    vert[2] = a;
    }

    compute_normals();
}

void Mesh :: apply_scale(){
    for(unsigned int i = 0; i < num_vertices; i++){
        vertices[i].x() = E[0] * vertices[i].x();
        vertices[i].y() = E[1] * vertices[i].y();
        vertices[i].z() = E[2] * vertices[i].z();
    }
    E[0] = E[1] = E[2] = 1.0f;
}

void Mesh :: save_to_off(const char *filename, bool flip){
    ofstream out(filename);

    out << "OFF" << endl;
    out << num_vertices << " " << num_triangles << " 0" << endl;
    for(unsigned int i = 0; i < num_vertices; i++){
        out << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << endl;
    }

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        if(!flip)
            out << "3 " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
        else
            out << "3 " << vert[0] << " " << vert[2] << " " << vert[1] << endl;
    }

    out.close();
}

void Mesh :: save_to_off_with_color(const char *filename){
    ofstream out(filename);

    out << "COFF" << endl;
    out << num_vertices << " " << num_triangles << " 0" << endl;
    for(unsigned int i = 0; i < num_vertices; i++){
        out << vertices[i].x()/E[0] << " " << vertices[i].y()/E[1] << " " << vertices[i].z()/E[2] << " " << vertices[i].red() << " " << vertices[i].green() << " " << vertices[i].blue() << endl;
    }

    for(unsigned int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        out << "3 " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
    }

    out.close();
}

float Mesh :: computeDiagonal(){
    float bbox[2][3] = { { 1.0E30F, 1.0E30F, 1.0E30F }, { -1.0E30F, -1.0E30F, -1.0E30F } };

    for (unsigned int i = 0; i < get_number_vertices(); i++) {
       if (vertices[i].x() < bbox[0][0]) bbox[0][0] = vertices[i].x();
       else if (vertices[i].x() > bbox[1][0]) bbox[1][0] = vertices[i].x();
       if (vertices[i].y() < bbox[0][1]) bbox[0][1] = vertices[i].y();
       else if (vertices[i].y() > bbox[1][1]) bbox[1][1] = vertices[i].y();
       if (vertices[i].z() < bbox[0][2]) bbox[0][2] = vertices[i].z();
       else if (vertices[i].z() > bbox[1][2]) bbox[1][2] = vertices[i].z();
    }

      // Setup initial viewing scale
      float dx = bbox[1][0] - bbox[0][0];
      float dy = bbox[1][1] - bbox[0][1];
      float dz = bbox[1][2] - bbox[0][2];

      return sqrt(dx*dx + dy*dy + dz*dz);
}

void Mesh :: decompose_connected_components(){
    vertex_components.clear();
    components.clear();

    for(unsigned int vertexSeed = 0; vertexSeed < num_vertices; ){
        //find first unselected vertex
        bool foundSeed = false;
        while(vertexSeed < num_vertices){
            if(!vertices[vertexSeed].is_selected()){
                foundSeed = true;
                break;
            }
            ++vertexSeed;
        }
        if(!foundSeed)
            break;

        components.resize(components.size() + 1);
        vector<Vertex*> activevertices;
        activevertices.push_back(&vertices[vertexSeed]);

        while(!activevertices.empty()){
            Vertex* v = activevertices.back();
            activevertices.pop_back();

            if(v->is_selected())
                continue;

            v->set_selected();
            components.back().push_back(v);
            set<unsigned int> av = v->get_vertices();

            for(set<unsigned int>::iterator it = av.begin();it!=av.end(); it++){
                if(!vertices[*it].is_selected())
                    activevertices.push_back(&vertices[*it]);
            }

        }

        ++vertexSeed;
    }

    vertex_components.resize(num_vertices, -1);

    for(unsigned int i = 0; i < components.size(); i++){
        for(unsigned int j = 0; j < components[i].size(); j++)
            vertex_components[components[i][j]->get_index()] = i;
    }

    //cout << "# components(decompose):" << components.size() << endl;
    marked_components.resize(components.size(), false);
}

Mesh* Mesh :: remove_marked_components(){
    set<unsigned int> vertReturned;
    set<unsigned int> faceReturned;

    for(unsigned int i = 0; i < components.size(); i++){
        if(!marked_components[i]){
            for(unsigned int j = 0; j < components[i].size(); j++){
                vertReturned.insert(components[i][j]->get_index());
                vector<unsigned int> faces = components[i][j]->get_faces();
                faceReturned.insert(faces.begin(), faces.end());
            }
        }
    }

    if(vertReturned.size()>0){
        Triangle* faces = triangles;

        set<unsigned int>::iterator it;
        for(it = faceReturned.begin(); it!=faceReturned.end(); it++){
            unsigned int faceInd = *it;
            unsigned int ind1 = faces[faceInd].get_vertex_at(0);
            unsigned int ind2 = faces[faceInd].get_vertex_at(1);
            unsigned int ind3 = faces[faceInd].get_vertex_at(2);

            vertReturned.insert(ind1);
            vertReturned.insert(ind2);
            vertReturned.insert(ind3);
        }

        Mesh* result = new Mesh();

        result->set_num_vertices(vertReturned.size());
        result->set_num_triangles(faceReturned.size());


        map<unsigned int, unsigned int> mapping;
        it = vertReturned.begin();

        unsigned int i = 0;
        while(it!=vertReturned.end()){
            unsigned int ind = *it;
            mapping.insert( make_pair(ind, i) );
            result->add_vertex(i, vertices[ind].x(), vertices[ind].y(), vertices[ind].z());
            it++;
            i++;
        }

        it = faceReturned.begin();
        i = 0;
        while(it!=faceReturned.end()){
            unsigned int faceInd = *it;
            unsigned int ind1 = faces[faceInd].get_vertex_at(0);
            unsigned int ind2 = faces[faceInd].get_vertex_at(1);
            unsigned int ind3 = faces[faceInd].get_vertex_at(2);

            result->add_triangle(i, mapping[ind1], mapping[ind2], mapping[ind3]);
            i++;
            it++;
        }

        result->compute_normals();
        return result;
     }
     return NULL;
}

void Mesh :: compute_fractured_components(){
    vertex_components.clear();
    components.clear();

    for(unsigned int vertexSeed = 0; vertexSeed < num_vertices; ){
        //find first unselected vertex
        bool foundSeed = false;
        while(vertexSeed < num_vertices){
            if(!vertices[vertexSeed].is_selected() && vertices[vertexSeed].getMarkFracture()!=0){
                foundSeed = true;
                break;
            }
            ++vertexSeed;
        }
        if(!foundSeed)
            break;

        components.resize(components.size() + 1);
        vector<Vertex*> activevertices;
        activevertices.push_back(&vertices[vertexSeed]);

        while(!activevertices.empty()){
            Vertex* v = activevertices.back();
            activevertices.pop_back();

            if(v->is_selected())
                continue;

            v->set_selected();
            components.back().push_back(v);
            set<unsigned int> av = v->get_vertices();

            for(set<unsigned int>::iterator it = av.begin();it!=av.end(); it++){
                if(!vertices[*it].is_selected() && vertices[*it].getMarkFracture()!=0)
                    activevertices.push_back(&vertices[*it]);
            }

        }

        if(components.back().size() < 200)
            components.pop_back();
        ++vertexSeed;
    }

    vertex_components.resize(num_vertices, -1);

    for(unsigned int i = 0; i < components.size(); i++){
        for(unsigned int j = 0; j < components[i].size(); j++)
            vertex_components[components[i][j]->get_index()] = i;
    }

    //cout << "# components:" << components.size() << endl;
}

void Mesh :: save_component(char* filename, unsigned int pos){
    set<unsigned int> vertReturned;
    set<unsigned int> faceReturned;

    for(unsigned int i = 0; i < components[pos].size(); i++){
        vertReturned.insert(components[pos][i]->get_index());
        vector<unsigned int> faces = components[pos][i]->get_faces();
        faceReturned.insert(faces.begin(), faces.end());
    }

    Triangle* faces = triangles;

    set<unsigned int>::iterator it;
    for(it = faceReturned.begin(); it!=faceReturned.end(); it++){
        unsigned int faceInd = *it;
        unsigned int ind1 = faces[faceInd].get_vertex_at(0);
        unsigned int ind2 = faces[faceInd].get_vertex_at(1);
        unsigned int ind3 = faces[faceInd].get_vertex_at(2);

        vertReturned.insert(ind1);
        vertReturned.insert(ind2);
        vertReturned.insert(ind3);
    }

    Mesh* result = new Mesh();

    result->set_num_vertices(vertReturned.size());
    result->set_num_triangles(faceReturned.size());


    map<unsigned int, unsigned int> mapping;
    it = vertReturned.begin();

    unsigned int i = 0;
    while(it!=vertReturned.end()){
        unsigned int ind = *it;
        mapping.insert( make_pair(ind, i) );
        result->add_vertex(i, vertices[ind].x(), vertices[ind].y(), vertices[ind].z());
        it++;
        i++;
    }

    it = faceReturned.begin();
    i = 0;
    while(it!=faceReturned.end()){
        unsigned int faceInd = *it;
        unsigned int ind1 = faces[faceInd].get_vertex_at(0);
        unsigned int ind2 = faces[faceInd].get_vertex_at(1);
        unsigned int ind3 = faces[faceInd].get_vertex_at(2);

        result->add_triangle(i, mapping[ind1], mapping[ind2], mapping[ind3]);
        i++;
        it++;
    }


    result->compute_normals();
    result->save_to_off(filename);
    delete result;
}

float* Mesh::constVertices() {
    p_Vertices = new float[num_vertices * 3];

    for(int i = 0; i < num_vertices; i++){
        p_Vertices[3 * i] = vertices[i].x();
        p_Vertices[3 * i + 1] = vertices[i].y();
        p_Vertices[3 * i + 2] = vertices[i].z();
    }

    return p_Vertices;
}

float* Mesh::constNormals() {
    p_Normals = new float[num_vertices * 3];

    for(int i = 0; i < num_vertices; i++){
        p_Normals[3 * i] = vertices[i].nx();
        p_Normals[3 * i + 1] = vertices[i].ny();
        p_Normals[3 * i + 2] = vertices[i].nz();
    }

    return p_Normals;
}

unsigned int* Mesh::constElements(){
    p_Elements = new unsigned int[num_triangles * 3];

    for(int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        p_Elements[3 * i] = vert[0];
        p_Elements[3 * i + 1] = vert[1];
        p_Elements[3 * i + 2] = vert[2];
    }

    return p_Elements;
}

}
