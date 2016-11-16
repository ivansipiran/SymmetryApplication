#include "object3d.h"

namespace SymArchData{

typedef Vertex Vector;

float determinant(Vector a, Vector b){
    return a.x()*(b.y() - b.z()) - a.y()*(b.x() - b.z()) + a.z()*(b.x() - b.y());
}

Object3D :: Object3D (){
    vertices = NULL;
    triangles = NULL;
    num_vertices = 0;
    num_triangles = 0;
    bounding_box = NULL;
    type_mesh = TRIANGLE_MESH;
    is_visible = true;
    E[0] = E[1] = E[2] = 1.0f;
}

Object3D :: ~Object3D(){
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
}

void Object3D :: compute_normals(){
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

void Object3D :: set_num_vertices(unsigned int nv){
    num_vertices = nv;
    vertices = new Vertex[num_vertices];
}

void Object3D :: set_num_triangles(unsigned int nt){
    num_triangles = nt;
    triangles = new Triangle[num_triangles];
}

void Object3D ::add_vertex(unsigned int pos, float x, float y, float z){
    vertices[pos].x() = x;
    vertices[pos].y() = y;
    vertices[pos].z() = z;
    vertices[pos].set_index(pos);
}

void Object3D :: set_vertex(unsigned int pos, float x, float y, float z){
    vertices[pos].x() = x;
    vertices[pos].y() = y;
    vertices[pos].z() = z;
}

void Object3D :: add_vertex(Vertex v, unsigned int pos)
{
    vertices[pos].x() = v.x();
    vertices[pos].y() = v.y();
    vertices[pos].z() = v.z();
    vertices[pos].set_index(pos);
    vertices[pos].set_weight(v.get_weight());
    vertices[pos].compute_color();
}

void Object3D :: add_triangle(unsigned int pos, unsigned int p1, unsigned int p2, unsigned int p3){
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

void Object3D :: add_normal(unsigned int pos, float nx, float ny, float nz){
    vertices[pos].set_normal(nx, ny, nz);
}

void Object3D :: add_color(unsigned int pos, unsigned int r, unsigned int g, unsigned int b){
    vertices[pos].set_color(r, g, b);
}

float Object3D::normalize_area(){
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
    //apply_scale();
    return sc;
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

} //namespace

