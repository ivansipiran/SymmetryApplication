#ifndef ANALYSISBAR_H
#define ANALYSISBAR_H

#include "SimpleMesh/mesh.h"

class AnalysisBar : public SimpleMesh::Mesh
{
protected:
    unsigned int mNumSectors;
    unsigned int mNumPartitions;
    float mRadius;

public:
    AnalysisBar();
    AnalysisBar(unsigned int pNumSectors, unsigned int pNumPartitions, float pRadius);
    ~AnalysisBar();

    void setNumberSectors(unsigned int pNumSectors){mNumSectors = pNumSectors;}
    void setNumPartitions(unsigned int pNumPartitions){mNumPartitions = pNumPartitions;}
    void setRadius(float pRadius){mRadius = pRadius;}

    void proceduralInit();
};

#endif // ANALYSISBAR_H
