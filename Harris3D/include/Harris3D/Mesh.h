#ifndef __MESH_H_
#define __MESH_H_

#include "harris_global.h"
#include <SimpleMesh/vertex.h>
#include "SimpleMesh/mesh.h"
#include <vector>
#include <cstring>

namespace Harris3D{

class HARRIS_API Mesh : public SimpleMesh::Mesh{

    bool* mark;
    bool* interest;
    unsigned int* depth;
    float* response;

  public:
      void cleanMesh();

      Mesh();
      ~Mesh();

      void initializeVertexInformation();
      void getNeighborhoodAtVertex(unsigned int index, int rad, std::vector<SimpleMesh::Vertex*>& V);
      int getRadiusAtVertex(unsigned int index, float radius, std::vector<SimpleMesh::Vertex*>& V);

      void processMaximumAtVertex(unsigned int index, int numRings);
      void setResponses(float* resp) {memcpy(response, resp, sizeof(float)*num_vertices);}
      bool getInterestAtVertex(unsigned int index){return interest[index];}

      float getResponseAt(unsigned int index){return response[index];}

      void deepCopy(SimpleMesh::Mesh* in);

};

}
#endif
