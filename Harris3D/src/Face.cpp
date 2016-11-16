#include "Harris3D/Face.h"

namespace Harris3D{

Face::Face(){

}

int Face::index(unsigned int v){
  if(vertices[0] == v) return 0;
  if(vertices[1] == v) return 1;
  if(vertices[2] == v) return 2;
  return -1;
}
}
