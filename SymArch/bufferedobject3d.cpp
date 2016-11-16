#include "bufferedobject3d.h"
#include "glwidget.h"
#include <QOpenGLShaderProgram>
#include <iostream>

using namespace std;

namespace SymArchData{

BufferedObject3D::BufferedObject3D()
    :m_vertexBuffer(QOpenGLBuffer::VertexBuffer),
      m_normalBuffer(QOpenGLBuffer::VertexBuffer),
      m_colorBuffer(QOpenGLBuffer::VertexBuffer),
      m_indexBuffer(QOpenGLBuffer::IndexBuffer)
{
    m_vertices = NULL;
    m_normals = NULL;
    m_indices = NULL;
    m_colors = NULL;
}

BufferedObject3D::~BufferedObject3D(){
    if(m_vertices){delete[] m_vertices; m_vertices = NULL;}
    if(m_normals){delete[] m_normals; m_normals = NULL;}
    if(m_indices){delete[] m_indices; m_indices = NULL;}
    if(m_colors){delete[] m_colors; m_colors = NULL;}

    if(m_vertexBuffer.isCreated())
        m_vertexBuffer.destroy();
    if(m_normalBuffer.isCreated())
        m_normalBuffer.destroy();
    if(m_colorBuffer.isCreated())
        m_colorBuffer.destroy();
    if(m_indexBuffer.isCreated())
        m_indexBuffer.destroy();
}

void BufferedObject3D::createBuffersFromObject(SimpleMesh::Mesh *mesh){

    glArea->makeCurrent();
    initializeOpenGLFunctions();
    num_vertices = mesh->get_number_vertices();
    SimpleMesh::Vertex* vertices = mesh->get_vertices();

    //Creating the vertex buffer
    m_vertices = new float[num_vertices * 3];
    for(int i = 0; i < num_vertices; i++){
        m_vertices[3 * i]       = vertices[i].x();
        m_vertices[3 * i + 1]   = vertices[i].y();
        m_vertices[3 * i + 2]   = vertices[i].z();
    }

    //Creating the normal buffer
    m_normals = new float[num_vertices * 3];
    for(int i = 0; i < num_vertices; i++){
        m_normals[3 * i]        = vertices[i].nx();
        m_normals[3 * i + 1]    = vertices[i].ny();
        m_normals[3 * i + 2]    = vertices[i].nz();
    }

    //Creating the index buffer
    num_triangles = mesh->get_number_triangles();
    SimpleMesh::Triangle* triangles = mesh->get_triangles();

    m_indices = new unsigned int[num_triangles * 3];
    for(int i = 0; i < num_triangles; i++){
        vector<unsigned int> vert = triangles[i].get_vertices();
        m_indices[3 * i]        = vert[0];
        m_indices[3 * i + 1]    = vert[1];
        m_indices[3 * i + 2]    = vert[2];
    }

    m_vao.create();
    QOpenGLVertexArrayObject::Binder vao_binder(&m_vao);

        m_vertexBuffer.create();
        m_vertexBuffer.bind();
        m_vertexBuffer.allocate(m_vertices, num_vertices * 3 * sizeof(GLfloat));

        m_normalBuffer.create();
        m_normalBuffer.bind();
        m_normalBuffer.allocate(m_normals, num_vertices * 3 * sizeof(GLfloat));

        m_indexBuffer.create();
        m_indexBuffer.bind();
        m_indexBuffer.allocate(m_indices, num_triangles * 3 * sizeof(GLuint));
        glArea->doneCurrent();
}

void BufferedObject3D::drawBufferedObject(QOpenGLShaderProgram* program){
    QOpenGLVertexArrayObject::Binder vao_binder(&m_vao);

    int material_color = program->uniformLocation("material_color");
    program->setUniformValue(material_color, color);

    m_vertexBuffer.bind();
        int vertexLocation = program->attributeLocation("v_coord");
        program->enableAttributeArray(vertexLocation);
        program->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(float) * 3);
    m_vertexBuffer.release();

    m_normalBuffer.bind();
        int normalLocation = program->attributeLocation("v_normal");
        program->enableAttributeArray(normalLocation);
        program->setAttributeBuffer(normalLocation, GL_FLOAT, 0, 3, sizeof(float) * 3);
    m_normalBuffer.release();

    m_indexBuffer.bind();
        glDrawElements(GL_TRIANGLES, num_triangles * 3, GL_UNSIGNED_INT, 0 );
    m_indexBuffer.release();

}
} //namespace
