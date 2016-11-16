#include <cassert>
#include "analysisbar.h"
#include <iostream>

AnalysisBar::AnalysisBar():
    Mesh()
{
    mNumSectors = 0;
    mNumPartitions = 0;
}

AnalysisBar::AnalysisBar(unsigned int pNumSectors, unsigned int pNumPartitions, float pRadius):
    Mesh()
{
    mNumSectors = pNumSectors;
    mNumPartitions = pNumPartitions;
    mRadius = pRadius;
}

AnalysisBar::~AnalysisBar(){

}

/*This method builds a cylinder based on the number of sectors and partitions we provide in the constructor
The pose of the cylinder is canonical, so it must be transformed after initialization*/
void AnalysisBar::proceduralInit(){
    //Number of vertices is numSectors * numPartitions + 2
    set_number_vertices(mNumSectors * mNumPartitions + mNumSectors + 2);

    //Number of triangles is numSectors x numPartitions x 2
    set_number_triangles(mNumSectors*mNumPartitions*2 + mNumSectors * 2);

    float angle = -(2*3.141592)/mNumSectors;
    float offset = 1.0/mNumPartitions;

    SimpleMesh::Vertex initPoint(mRadius, 0.0, 0.0);
    //We start with point (1,0,0) and start the population transforming this point
    int currentVertex = 0;
    add_vertex(currentVertex++, initPoint);

    add_color(currentVertex-1, 1.0, 0.0, 0.0);

    //Populate the first level
    SimpleMesh::Vertex auxPoint;
    for(int i = 0; i < (mNumSectors - 1); i++){
        auxPoint.setX(initPoint.x()*cos(angle*(i+1)) - initPoint.y()*sin(angle*(i+1)));
        auxPoint.setY(initPoint.x()*sin(angle*(i+1)) + initPoint.y()*cos(angle*(i+1)));
        auxPoint.setZ(initPoint.z());

        add_vertex(currentVertex++, auxPoint);
        add_color(currentVertex-1, 1.0, 0.0, 0.0);
    }

    //Populate the next levels
    for(int i = 0; i < mNumPartitions; i++){
        for(int j = 0; j < mNumSectors; j++){
            auxPoint.setX(vertices[j].x());
            auxPoint.setY(vertices[j].y());
            auxPoint.setZ(vertices[j].z() + offset*(i+1));

            add_vertex(currentVertex++, auxPoint);
            add_color(currentVertex-1, 1.0, 0.0, 0.0);
        }
    }


    //The last two points are the center of cylinder's circles
    add_vertex(currentVertex++, 0.0, 0.0, 0.0);
    add_color(currentVertex-1, 1.0, 0.0, 0.0);

    add_vertex(currentVertex++, 0.0, 0.0, 1.0);
    add_color(currentVertex-1, 1.0, 0.0, 0.0);

    std::cout << currentVertex  << " -> " << num_vertices << endl;
    //Check number of vertices
    assert(currentVertex == num_vertices);

    int currentTriangle = 0;
    //Now it is the turn for the triangle population
    for(int i = 0; i < mNumPartitions; i++){
        for(int j = 0; j < (mNumSectors - 1); j++){
            int auxVal = i * mNumSectors + j;
            add_triangle(currentTriangle++, auxVal, auxVal + mNumSectors + 1, auxVal + 1);
            add_triangle(currentTriangle++, auxVal, auxVal + mNumSectors, auxVal + mNumSectors + 1);
        }
        int auxVal = (i + 1) * mNumSectors;
        add_triangle(currentTriangle++, auxVal - 1, auxVal, i * mNumSectors);
        add_triangle(currentTriangle++, auxVal - 1, auxVal - 1 + mNumSectors, auxVal);
    }

    //Last triangles are the top and bottom faces
    for(int i = 0; i <  mNumSectors; i++){
        add_triangle(currentTriangle++, currentVertex - 2, i, (i + 1)%mNumSectors);
        add_triangle(currentTriangle++, currentVertex - 1, mNumSectors * mNumPartitions + (i + 1)%mNumSectors, mNumSectors * mNumPartitions + i);
    }

    compute_normals();
}
