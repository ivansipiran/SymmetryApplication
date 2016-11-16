#ifndef __VERTEX_H_
#define __VERTEX_H_

#include <set>
#include <vector>
#include <cmath>
using namespace std;

namespace SymArchData{
/*!
  It stores vertex information. It contains geometric and topological information.
*/
class Vertex{
    private:
        float                       v[3];                /**< Coordinates of a 3D vertex*/
        float                       normal[3];           /**< Normal of the 3D vertex*/
        vector<unsigned int>        faces;               /**< Indices of incident faces*/
        set<unsigned int>           adjacentVertices;    /**< Set of adjacent vertices*/
        unsigned int                index;               /**< Unique identifier for a vertex*/

        float          weight;              /**< Value that stores additional information for the vertex*/
        unsigned int    r;                   /**< Red component of the vertex*/
        unsigned int    g;                   /**< Green component of the vertex*/
        unsigned int    b;                   /**< Blue component of the vertex*/
        float           rgb;                 /**< Real number that represent the color in one variable*/
        bool            isInterest;          /**< Flag that indicates the interest value of the vertex*/

        float          response;            /**< Value for storing the feature response of the vertex*/

        unsigned int mark_fracture;
        unsigned int mark_intact;
        unsigned int confidence;

    public:

       //! Base constructor
       /*!
         Initialize the basic fields of the object
        */
       Vertex(){v[0] = v[1] = v[3] = 0.0; isInterest=false;}

       //! Constructor
       /*!
         Constructor to initialize the coordinates of the vertex
         \param x1 Coordinate x
         \param y1 Coordinate y
         \param z1 Coordinate z
        */
       Vertex(float x1, float y1, float z1) {v[0] = x1; v[1] = y1; v[2] = z1;isInterest=false;}

       //! Access function
       /*!
          Get the x coordinate
        */
       float x() const {return v[0];}

       //! Access function
       /*!
          Get the y coordinate
        */
       float y() const {return v[1];}

       //! Access function
       /*!
          Get the z coordinate
        */
       float z() const {return v[2];}


       //! Access function with reference
       /*!
          Get the x coordinate (reference)
        */
       float& x() {return v[0];}

       //! Access function with reference
       /*!
          Get the y coordinate (reference)
        */
       float& y() {return v[1];}

       //! Access function with reference
       /*!
          Get the z coordinate (reference)
        */
       float& z() {return v[2];}

       //! Assignment function
       /*!
          \param x X-coordinate
        */
       void setX(float x){v[0] = x;}

       //! Assignment function
       /*!
          \param y Y-coordinate
        */
       void setY(float y){v[1] = y;}

       //! Assignment function
       /*!
          \param z Z-coordinate
        */
       void setZ(float z){v[2] = z;}

       //! Assignment function for normals
       /*!
          \param nx Component x of the normal
          \param ny Component y of the normal
          \param nz Component z of the normal
        */
       void set_normal(float nx, float ny, float nz) {normal[0] = nx; normal[1] = ny; normal[2] = nz;}

       //! Access function
       /*!
          Get the component x of the normal
        */
       float nx() const {return normal[0];}

       //! Access function
       /*!
          Get the component y of the normal
        */
       float ny() const {return normal[1];}

       //! Access function
       /*!
          Get the component z of the normal
        */
       float nz() const {return normal[2];}

       //! Access function with reference
       /*!
          Get the component x of normal (reference)
        */
       float& nx() {return normal[0];}

       //! Access function with reference
       /*!
          Get the component y of normal (reference)
        */
       float& ny() {return normal[1];}

       //! Access function with reference
       /*!
          Get the component < of normal (reference)
        */
       float& nz() {return normal[2];}

       //! Access function for incident faces
       /*!
         Get the indicent faces for a vertex
        */
       vector<unsigned int>& get_faces() {return faces;}
       set<unsigned int>& get_vertices(){return adjacentVertices;}

       //! Assignment function for incident faces
       /*!
         Associate the vertex with a face
         \param f Face's index to be associated with the vertex
        */
       void add_face(unsigned int f) {faces.push_back(f);}

       //! Assignment function for index
       /*!
         Associate an index with the vertex
         \param i Index to be associated with the vertex
        */
       void set_index(unsigned int i){index = i;}

       //! Access function for the index
       /*!
         Get the index of a vertex
        */
       unsigned int get_index(){return index;}

       //! Assignment function for weight
       /*!
         Associate a weight with the vertex
         \param w Weight to be associated with the vertex
        */
       void set_weight(float w) {weight = w;}

       //! Access function for the weight
       /*!
         Get the weight of a vertex
        */
       float get_weight(){ return weight;}

       //! Access function for the float rgb
       /*!
         Get the float color of the vertex
        */
       float get_rgb() {return rgb;}

       //! Access function for the rgb component
       /*!
         Get the red component of the vertex
        */
       unsigned int red() {return r;}

       //! Access function for the rgb component
       /*!
         Get the green component of the vertex
        */
        unsigned int green() {return g;}

        //! Access function for the rgb component
        /*!
          Get the  blue component of the vertex
         */
        unsigned int blue() {return b;}

       //! Compute a RGB color depending on the weight value
       void compute_color();

       //! Assignment function for response
       /*!
         Associate a response with the vertex
         \param resp Response to be associated with the vertex
        */
       void setResponse(float resp){response = resp;}

       //! Access function for the response
       /*!
         Get the  response of the vertex
        */
       float getResponse(){return response;}

       //! Computes the 1-neighborhood local maxima using the response field
       /*!
         This function compares the response of the vertex with the responses of the neighbor vertices. If the vertex's response is greater
         than the neighbor's responses, the mark of interest of the vertex is activated.
         \param vertices The set of vertices
        */
       void processMaximum(Vertex* vertices);

       //! Access function for the mark of interest
       /*!
         Get the  mark of interest of the vertex
        */
       bool getInterest(){return isInterest;}
       void setInterest(bool f){ isInterest = f;}

       //! Assignment function for adjacent vertices
       /*!
         Associate the vertex with its adjacent vertex
         \param vertex Vertex's index to be associated with this vertex
        */
       void addVertex(unsigned int vertex){ adjacentVertices.insert(vertex);}

       float dot(Vertex v1){ return (this->x()*v1.x() + this->y()*v1.y() + this->z()*v1.z());}

       void setMarkFracture(unsigned int mark){mark_fracture = mark;}
       void setMarkIntact(unsigned int mark){mark_intact = mark;}
       void setMarkConfidence(unsigned int conf){confidence = conf;}

       unsigned int getMarkFracture(){return r;}
       unsigned int getMarkIntact(){return b;}
       unsigned int getMarkConfidence(){return g;}

       void set_color(unsigned int c1, unsigned int c2, unsigned int c3){r = c1; g = c2; b = c3;}

       bool is_selected(){return isInterest;}
       void set_selected(){isInterest = true;}

       void normalize(){
           float mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
           v[0] = v[0]/mag;
           v[1] = v[1]/mag;
           v[2] = v[2]/mag;

       }

};

typedef Vertex Point3D;
typedef Vertex Vector3D;
typedef Vertex Vector;

}
#endif
