#include <cstdint>
#include "Vertex.h"

namespace SymArchData{

typedef struct {
    float r,g,b;
} COLOUR;

void Vertex :: compute_color(){
    float vmin = 0.0;
    float vmax = 1.0;
    float v = weight;

    COLOUR c = {1.0,1.0,1.0}; // white
    float dv;

    if (v < vmin)
       v = vmin;
    if (v > vmax)
       v = vmax;
    dv = vmax - vmin;

    if (v < (vmin + 0.25 * dv)) {
       c.r = 0;
       c.g = 4 * (v - vmin) / dv;
    } else if (v < (vmin + 0.5 * dv)) {
       c.r = 0;
       c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
    } else if (v < (vmin + 0.75 * dv)) {
       c.r = 4 * (v - vmin - 0.5 * dv) / dv;
       c.b = 0;
    } else {
       c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
       c.b = 0;
    }

    r = (unsigned int)(255*c.r);
    g = (unsigned int)(255*c.g);
    b = (unsigned int)(255*c.b);

    uint32_t urgb;
    urgb = r<<16 | g<<8 | b;
    rgb=*(reinterpret_cast<float *>(&urgb));
}

void Vertex::processMaximum(Vertex* vertices){
        set<unsigned int> :: iterator it;
        for(it = adjacentVertices.begin(); it!=adjacentVertices.end(); it++){
            Vertex* v1 = &vertices[*it];
            if(v1!=this){
                if(response < v1->getResponse())
                    return;
            }
        }
    isInterest = true;
}
}
