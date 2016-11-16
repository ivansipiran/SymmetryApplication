#ifndef __VERTEX_H_
#define __VERTEX_H_

#include "harris_global.h"
#include "Face.h"
#include "SimpleMesh/vertex.h"
#include <iostream>
#include <vector>
#include <queue>

#include <CGAL/basic.h>
#include <CGAL/Search_traits.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace Harris3D{

class Face;


class HARRIS_API Vertex : public SimpleMesh::Vertex {
 private:

      bool mark;
      bool isInterest;
      int depth;
      double response;
      int numNeighbors;



 public:

      Vertex() {coord[0] = coord[1] = coord[2] = normal[0] = normal[1] = normal[2] = 0.0; mark = false; depth = 0;isInterest = 0;}
      Vertex(float x1, float y1, float z1) {coord[0] = x1; coord[1] = y1; coord[2] = z1; mark = false; depth = 0;isInterest = 0;}

      inline bool operator==(const Vertex& p) const { return (x() == p.x()) && (y() == p.y()) && (z() == p.z()); }
      inline bool  operator!=(const Vertex& p) const { return ! (*this == p); }

      friend std::ostream& operator<<(std::ostream& out, Vertex& point);

      inline bool isMarked(){return mark;}
      inline void setMark(bool mark1){mark = mark1;}
      inline int getDepth(){ return depth;}
      inline void setDepth(int dep){ depth = dep;}
      inline double getResponse(){return response;}
      inline void setResponse(double resp){response = resp;}
      inline bool getInterest(){return isInterest;}

      void getNeighborhood(int rad, std::vector<Vertex*>& V, Vertex* vertices);
      int getRadius(Vertex* vertices, float radius, std::vector<Vertex*>& V);

      void processMaximum(Vertex* vertices, int numRings);

      void getPatch(Vertex* vertices, std::vector<unsigned int> indices, std::set<unsigned int>& returned, std::set<unsigned int>& faceR, float radius, Vertex center);


      const float* get_coord(){return coord;}
	float distanceL2(Vertex* v1);
  };
}

namespace CGAL {

  template <>
  struct Kernel_traits<Harris3D::Vertex> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}


/*struct Construct_coord_iterator {
  const float* operator()(const Harris3D::Vertex& p) const
  { return static_cast<const float*>(p.get_coord()); }

  const float* operator()(const Harris3D::Vertex& p, int)  const
  { return static_cast<const float*>(p.get_coord()+3); }
};*/
#endif
