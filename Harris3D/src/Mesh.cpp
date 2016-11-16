#include <map>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <cassert>
#include <queue>
#include <cfloat>
#include "Harris3D/Mesh.h"
#include "Util/util.h"
//#include "Util/Clock.h"

namespace Harris3D{

Mesh::Mesh():
    SimpleMesh::Mesh(){
       diag = -1.0;
       mark = 0;
       interest = 0;
       depth = 0;
       response = 0;
}


void Mesh :: initializeVertexInformation(){
    mark = new bool[num_vertices];
    interest = new bool[num_vertices];
    depth = new unsigned int[num_vertices];
    response = new float[num_vertices];

    memset(mark, 0, sizeof(bool)*num_vertices);
    memset(interest, 0, sizeof(bool)*num_vertices);
    memset(depth, 0, sizeof(unsigned int)*num_vertices);
    memset(response, 0, sizeof(float)*num_vertices);
}

int Mesh :: getRadiusAtVertex(unsigned int index, float radius, std::vector<SimpleMesh::Vertex *> &V){
    std::vector<SimpleMesh::Vertex*> marked; //Store the marked vertices
    std::map<unsigned int, float> distances; //Store the distances relatives to the current vertex
    std::map<unsigned int, SimpleMesh::Vertex*> markedRing; //Elements in a new ring
    std::queue<SimpleMesh::Vertex*> Q;
    float maxDistance = 0.0;
    int rad = -1;

    Q.push(&vertices[index]);
    mark[index] = true;
    markedRing.insert(std::pair<unsigned int, SimpleMesh::Vertex*>(index, &vertices[index]));

    distances[index] = 0.0;

    while(!Q.empty()){
        SimpleMesh::Vertex* v0 = Q.front();
        Q.pop();

        int dep = depth[v0->get_index()];
        if(dep != rad){ //First vertex in the new ring
            std::map<unsigned int, SimpleMesh::Vertex*>::iterator it;
            float max = 0.0;

            //Mark the previous ring
            for(it = markedRing.begin(); it!=markedRing.end(); it++){
                SimpleMesh::Vertex* mar = (*it).second;
                //mar->setMark(true);
                mark[mar->get_index()] = true;
                marked.push_back(mar);
                V.push_back(new SimpleMesh::Vertex(mar->x(), mar->y(), mar->z()));
                if(distances[(*it).first] > max)
                    max = distances[(*it).first];
            }

            rad++;
            markedRing.clear();
            maxDistance = max;
            if(maxDistance > radius)
                break;
        }

        std::set<unsigned int> listVertices = v0->get_vertices();
        std::set<unsigned int> :: iterator it;

        for(it = listVertices.begin(); it!=listVertices.end(); it++){
                SimpleMesh::Vertex* v1 = &vertices[*it];
                if(!mark[v1->get_index()]){
                    if(distances[v1->get_index()] == 0.0){ //Distance is not set
                        Q.push(v1);
                        //v1->setDepth(dep + 1);
                        depth[v1->get_index()] = dep + 1;
                    }
                    markedRing.insert(std::pair<unsigned int, SimpleMesh::Vertex*>(v1->get_index(), v1));
                    float dist = v0->distanceL2(v1);
                    float newDistance = distances[v0->get_index()] + dist;
                    if(distances[v1->get_index()] == 0.0){ //First time on this vertex
                        distances[v1->get_index()] = newDistance;
                    }else if(newDistance  < distances[v1->get_index()]){
                        distances[v1->get_index()] = newDistance;
                    }
                }
            }
        //}
    }

    if(!markedRing.empty()){
            std::map<unsigned int, SimpleMesh::Vertex*>::iterator it;
            float max = 0.0;


            for(it = markedRing.begin(); it!=markedRing.end(); it++){
                SimpleMesh::Vertex* mar = (*it).second;
                //mar->setMark(true);
                mark[mar->get_index()] = true;
                marked.push_back(mar);
                V.push_back(new SimpleMesh::Vertex(mar->x(), mar->y(), mar->z()));
                if(distances[(*it).first] > max)
                    max = distances[(*it).first];
            }

            rad++;
            markedRing.clear();
            maxDistance = max;

    }

    //Unmark all vertices
    //std::vector<SimpleMesh::Vertex*>::iterator ini = marked.begin();
    //while(ini  < marked.end()){
    //	(*ini)->setMark(false);
    //	(*ini)->setDepth(0);
    //	ini++;
    //}
    memset(mark, 0, sizeof(bool)*num_vertices);
    memset(depth, 0, sizeof(unsigned int)*num_vertices);

    return rad;
}

void Mesh::getNeighborhoodAtVertex(unsigned int index, int rad, std::vector<SimpleMesh::Vertex*>& V){
    std::queue<SimpleMesh::Vertex*> Q;
    std::vector<SimpleMesh::Vertex*> marked;

    //Q.push(this);
    Q.push(&vertices[index]);
    //this->setMark(true);
    //this->setDepth(0);
    mark[index] = true;
    depth[index] = 0;
    marked.push_back(&vertices[index]);

    while(!Q.empty()){
        SimpleMesh::Vertex* v0 = Q.front();
        Q.pop();
        V.push_back(new SimpleMesh::Vertex(v0->x(), v0->y(), v0->z())); //Indeed, copy vertex information rather than return the same vertex

        //int dep = v0->getDepth();
        int dep = depth[v0->get_index()];
        if(dep <= rad){
            std::set<unsigned int> listVertices = v0->get_vertices();
            std::set<unsigned int> :: iterator it;
            for(it = listVertices.begin(); it!=listVertices.end(); it++){
                SimpleMesh::Vertex* v1 = &vertices[*it];
                    //if(!v1->isMarked()){
                    if(!mark[v1->get_index()]){
                        Q.push(v1);
                        //v1->setMark(true);
                        //v1->setDepth(dep + 1);
                        mark[v1->get_index()] = true;
                        depth[v1->get_index()] = dep + 1;
                        marked.push_back(v1);
                    }
                }
            }
    }

    //std::vector<Vertex*>::iterator ini = marked.begin();

    //while(ini<marked.end()){
    //	(*ini)->setMark(false);
    //	(*ini)->setDepth(0);
    //	ini++;
    //}
    memset(mark, 0, sizeof(bool)*num_vertices);
    memset(depth, 0, sizeof(unsigned int)*num_vertices);

}

void Mesh::processMaximumAtVertex(unsigned int index, int numRings){
        std::set<unsigned int> :: iterator it;
        for(it = vertices[index].get_vertices().begin(); it!=vertices[index].get_vertices().end(); it++){
            SimpleMesh::Vertex* v1 = &vertices[*it];
            if(v1->get_index()!=index){
                if(response[index] < response[v1->get_index()])
                    return;
            }
        }
    interest[index] = true;
}

void Mesh::cleanMesh(){
    if(vertices){
        delete[] vertices;
        vertices = 0;
        num_vertices = 0;
    }

    if(triangles){
        delete[] triangles;
        triangles = 0;
        num_triangles = 0;
    }

    if(mark){
        delete[] mark;
        mark = 0;
    }
    if(interest){
        delete[] interest;
        interest = 0;
    }
    if(depth){
        delete[] depth;
        depth = 0;
    }
    if(response){
        delete[] response;
        response = 0;
    }
    diag = -1.0;
}

Mesh::~Mesh(){
    	cleanMesh();
}

void Mesh::deepCopy(SimpleMesh::Mesh *in){

    //Inicializar vertices y triangulos
    this->set_number_vertices(in->get_number_vertices());
    this->set_number_triangles(in->get_number_triangles());

    //Copiar vertices
    for(unsigned int i = 0; i < in->get_number_vertices(); i++){
        this->add_vertex(i, in->get_vertices()[i].x(), in->get_vertices()[i].y(), in->get_vertices()[i].z());
        if(in->has_normals())
            this->add_normal(i, in->get_vertices()[i].nx(), in->get_vertices()[i].ny(), in->get_vertices()[i].nz());
        if(in->has_color())
            this->add_color(i, in->get_vertices()[i].red(), in->get_vertices()[i].green(), in->get_vertices()[i].blue());
    }

    //Copiar triangulos
    for(unsigned int i = 0; i < in->get_number_triangles(); i++){
        this->add_triangle(i, in->get_triangles()[i].get_vertex_at(0), in->get_triangles()[i].get_vertex_at(1), in->get_triangles()[i].get_vertex_at(2));
    }

    //Copiar Flags
    if(in->has_normals())
        this->flag_normals = true;
    if(in->has_color())
        this->flag_color = true;

}
}


