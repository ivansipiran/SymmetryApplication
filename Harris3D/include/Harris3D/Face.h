#ifndef __FACE_H_
#define __FACE_H_

#include "harris_global.h"
#include <SimpleMesh/triangle.h>
#include <vector>

namespace Harris3D{

class HARRIS_API Face : public SimpleMesh::Triangle{

 public:
     Face();
     int index(unsigned int v);
};

}
#endif
