#ifndef OBJECTCONTAINER_H
#define OBJECTCONTAINER_H

#include <QVector3D>
#include <QOpenGLShaderProgram>
#include "SimpleMesh/mesh.h"
#include "bufferedobject3d.h"

//This class stores the objects in the application. This object contains two collections, the original objects and
//the corresponding buffered objects
class ObjectContainer
{
    vector<SimpleMesh::Mesh*> m_objects;
    vector<SymArchData::BufferedObject3D*> m_bufferedObjects;

    float bbox[2][3];
    float scale;
    QVector3D center;

public:
    ObjectContainer();
    ~ObjectContainer();

    void addObject(SimpleMesh::Mesh* new_mesh, SymArchData::BufferedObject3D* new_buffered_object);
    void drawObjects(QOpenGLShaderProgram* program);
    float getScale(){return scale;}
    QVector3D getCenter(){return center;}
    void removeObject(int index);
    void recomputeScaleAndCenter();

    void drawBoundingBox(QOpenGLShaderProgram* program);
    void clearData();
    int getNumberObjects(){return m_objects.size();}
    SimpleMesh::Mesh* getObject(int index){return m_objects[index];}
};

#endif // OBJECTCONTAINER_H
