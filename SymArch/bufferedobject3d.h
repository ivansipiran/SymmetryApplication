#ifndef BUFFEREDOBJECT3D_H
#define BUFFEREDOBJECT3D_H

#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QVector3D>
#include "SimpleMesh/mesh.h"

/**
 * The class BufferedObject3D stores complementary information for a mesh
 *
 * Basically, it stores the information related to linear buffers and opengl stuff
 */
#include <QOpenGLFunctions>

QT_FORWARD_DECLARE_CLASS(GLWidget)
QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

namespace SymArchData{

class BufferedObject3D : protected QOpenGLFunctions
{
private:
    float*          m_vertices;
    float*          m_normals;
    float*          m_colors;
    unsigned int*   m_indices;

    QOpenGLBuffer   m_vertexBuffer;
    QOpenGLBuffer   m_normalBuffer;
    QOpenGLBuffer   m_colorBuffer;
    QOpenGLBuffer   m_indexBuffer;

    QOpenGLVertexArrayObject m_vao;
    int num_vertices;
    int num_triangles;
    GLWidget* glArea;
    QVector3D color;

public:
    BufferedObject3D();
    ~BufferedObject3D();

    //This function transforms the data from a storage Mesh into a buffered object
    //It also initializes the VAO and VBO's
    void createBuffersFromObject(SimpleMesh::Mesh* mesh);
    void drawBufferedObject(QOpenGLShaderProgram* program);
    void setGLWidget(GLWidget* glWidget){glArea = glWidget;}
    void setColor(QVector3D v){color.setX(v.x()); color.setY(v.y()); color.setZ(v.z());}
};

}
#endif // BUFFEREDOBJECT3D_H
