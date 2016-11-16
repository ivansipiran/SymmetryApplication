#ifndef VERTEX_H
#define VERTEX_H

#include "SimpleMesh/simplemesh_global.h"
#include <vector>
#include <set>
#include <cmath>

using namespace std;

namespace SimpleMesh{

//! A class to store vertex information
/*!
 * Vertex stores the basic information for a 3D vertex in a mesh.
*/

class SIMPLEMESH_API Vertex
{

protected:
    float                   coord[3];           //!< 3D coordinates
    float                   normal[3];          //!< Vertex normal
    float                   color[3];           //!< Vertex color
    vector<unsigned int>    triangles;          //!< Indices of incident triangles
    set<unsigned int>       adjacentVertices;   //!< Indices of adjacent vertices
    unsigned int            index;              //!< Index of the vertex (a unique identifier)

    //Store information on barycentric representation
    unsigned int                     numTriangle;        //!< Index of incident triangle (any of them)
    unsigned int                     posBarycentric;     //!< Barycentric position

public:
    //! Default constructor
    /*! This constructor sets the coord and normal to 0.0
    */
    Vertex(){coord[0] = coord[1] = coord[2] = normal[0] = normal[1] = normal[2] = 0.0; color[0]=color[1]=color[2]=0.0;}

    //! Coord constructor
    /*! This constructor sets the 3D coordinates of a vertex
    */
    Vertex(float x, float y, float z){coord[0] = x; coord[1] = y; coord[2] = z;}

    //! Const getter for the x-coordinate
    float x() const {return coord[0];}

    //! Const getter for the y-coordinate
    float y() const {return coord[1];}

    //! Const getter for the z-coordinate
    float z() const {return coord[2];}

    //! Ref getter for the x-coordinate
    float& x() {return coord[0];}

    //! Ref getter for the y-coordinate
    float& y() {return coord[1];}

    //! Ref getter for the z-coordinate
    float& z() {return coord[2];}

    //! Setter for the x-coordinate
    void setX(float x){coord[0] = x;}

    //! Setter for the y-coordinate
    void setY(float y){coord[1] = y;}

    //! Setter for the z-coordinate
    void setZ(float z){coord[2] = z;}

    //! Setter for the normal
    void set_normal(float nx, float ny, float nz){normal[0] = nx; normal[1] = ny; normal[2] = nz;}

    //! Const getter for the normal (X coord)
    float nx() const {return normal[0];}

    //! Const getter for the normal (Y coord)
    float ny() const {return normal[1];}

    //! Const getter for the normal (Z coord)
    float nz() const {return normal[2];}

    //! Ref getter for the normal (X coord)
    float& nx() {return normal[0];}

    //! Ref getter for the normal (Y coord)
    float& ny() {return normal[1];}

    //! Ref getter for the normal (Z coord)
    float& nz() {return normal[2];}

    //! Setter for the color
    void set_color(float cx, float cy, float cz){color[0] = cx; color[1] = cy; color[2] = cz;}

    //! Const getter for the color (red)
    float red() const {return color[0];}

    //! Const getter for the color (green)
    float green() const {return color[1];}

    //! Const getter for the color (blue)
    float blue() const {return color[2];}

    //! Ref getter for the color (red)
    float& red() {return color[0];}

    //! Ref getter for the color (green)
    float& green() {return color[1];}

    //! Ref getter for the color (blue)
    float& blue() {return color[2];}


    //! Ref getter of the incident triangles
    vector<unsigned int>& get_triangles() {return triangles;}

    //! Ref getter of the adjacent vertices
    set<unsigned int>& get_vertices() {return adjacentVertices;}

    //! Adds an incident triangle to the current vertex
    void add_triangle(unsigned int tri) {triangles.push_back(tri);}

    //! Setter of the vertex's index
    void set_index(unsigned int ind) {index = ind;}

    //! Getter of the vertex's index
    unsigned int get_index() {return index;}

    //! Adds an adjacent vertex to the current one
    void add_vertex(unsigned int vert) {adjacentVertices.insert(vert);}

    void set_vertex(float x, float y, float z){ coord[0] = x; coord[1] = y; coord[2] = z;}

    inline void set_num_triangle(unsigned int num){numTriangle = num;}
    inline unsigned int get_num_triangle(){return numTriangle;}

    inline void set_pos_barycentric(int pos){posBarycentric = pos;}
    inline unsigned int get_pos_barycentric(){return posBarycentric;}

    float distanceL2(Vertex* v1){
        float val = (x()-v1->x())*(x()-v1->x()) + (y()-v1->y())*(y()-v1->y()) + (z()-v1->z())*(z()-v1->z());
        return sqrt(val);
    }

};

} //namespace SimpleMesh

#endif // VERTEX_H
