#include "objectcontainer.h"
#include "SimpleMesh/vertex.h"
#include <iostream>

ObjectContainer::ObjectContainer()
{
    bbox[0][0] = bbox[0][1] = bbox[0][2] = 1.0E30F;
    bbox[1][0] = bbox[1][1] = bbox[1][2] = -1.0E30F;

}

ObjectContainer::~ObjectContainer(){
    clearData();
}

void ObjectContainer::clearData(){
    for(int i = 0; i < m_objects.size(); i++){
        delete m_objects[i];
        delete m_bufferedObjects[i];
    }
    m_objects.clear();
    m_bufferedObjects.clear();
    bbox[0][0] = bbox[0][1] = bbox[0][2] = 1.0E30F;
    bbox[1][0] = bbox[1][1] = bbox[1][2] = -1.0E30F;
}

void ObjectContainer::removeObject(int index){
    SimpleMesh::Mesh* auxMesh = m_objects[index];
    SymArchData::BufferedObject3D* auxBO = m_bufferedObjects[index];
    delete auxMesh;
    delete auxBO;

    vector<SimpleMesh::Mesh*>::iterator itObject = m_objects.begin();
    vector<SymArchData::BufferedObject3D*>::iterator itBO = m_bufferedObjects.begin();

    m_objects.erase(itObject + index);
    m_bufferedObjects.erase(itBO + index);
}

void ObjectContainer::recomputeScaleAndCenter(){
    bbox[0][0] = bbox[0][1] = bbox[0][2] = 1.0E30F;
    bbox[1][0] = bbox[1][1] = bbox[1][2] = -1.0E30F;

    for(int j = 0; j < m_objects.size(); j++){
        //Update the bounding box
        SimpleMesh::Vertex* verts = m_objects[j]->get_vertices();
        for (int i = 0; i < m_objects[j]->get_number_vertices(); i++) {
            if (verts[i].x() < bbox[0][0]) bbox[0][0] = verts[i].x();
            else if (verts[i].x() > bbox[1][0]) bbox[1][0] = verts[i].x();
            if (verts[i].y() < bbox[0][1]) bbox[0][1] = verts[i].y();
            else if (verts[i].y() > bbox[1][1]) bbox[1][1] = verts[i].y();
            if (verts[i].z() < bbox[0][2]) bbox[0][2] = verts[i].z();
            else if (verts[i].z() > bbox[1][2]) bbox[1][2] = verts[i].z();
        }
    }

    // Setup initial viewing scale
    float dx = bbox[1][0] - bbox[0][0];
    float dy = bbox[1][1] - bbox[0][1];
    float dz = bbox[1][2] - bbox[0][2];
    scale = 2.0 / sqrt(dx*dx + dy*dy + dz*dz);

    center[0] = 0.5 * (bbox[1][0] + bbox[0][0]);
    center[1] = 0.5 * (bbox[1][1] + bbox[0][1]);
    center[2] = 0.5 * (bbox[1][2] + bbox[0][2]);

}

void ObjectContainer::addObject(SimpleMesh::Mesh *new_mesh, SymArchData::BufferedObject3D *new_buffered_object)
{
    m_objects.push_back(new_mesh);
    m_bufferedObjects.push_back(new_buffered_object);

    //Update the bounding box
    SimpleMesh::Vertex* verts = new_mesh->get_vertices();
    for (int i = 0; i < new_mesh->get_number_vertices(); i++) {
        if (verts[i].x() < bbox[0][0]) bbox[0][0] = verts[i].x();
        else if (verts[i].x() > bbox[1][0]) bbox[1][0] = verts[i].x();
        if (verts[i].y() < bbox[0][1]) bbox[0][1] = verts[i].y();
        else if (verts[i].y() > bbox[1][1]) bbox[1][1] = verts[i].y();
        if (verts[i].z() < bbox[0][2]) bbox[0][2] = verts[i].z();
        else if (verts[i].z() > bbox[1][2]) bbox[1][2] = verts[i].z();
    }

    // Setup initial viewing scale
    float dx = bbox[1][0] - bbox[0][0];
    float dy = bbox[1][1] - bbox[0][1];
    float dz = bbox[1][2] - bbox[0][2];
    scale = 2.0 / sqrt(dx*dx + dy*dy + dz*dz);

    center[0] = 0.5 * (bbox[1][0] + bbox[0][0]);
    center[1] = 0.5 * (bbox[1][1] + bbox[0][1]);
    center[2] = 0.5 * (bbox[1][2] + bbox[0][2]);

}

void ObjectContainer::drawObjects(QOpenGLShaderProgram* program){
    for(int i = 0; i < m_bufferedObjects.size(); i++)
        m_bufferedObjects[i]->drawBufferedObject(program);
}

void ObjectContainer::drawBoundingBox(QOpenGLShaderProgram *program){
    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();

    GLfloat vertices[]={
        bbox[0][0], bbox[0][1], bbox[0][2], 1.0,
        bbox[1][0], bbox[0][1], bbox[0][2], 1.0,
        bbox[1][0], bbox[1][1], bbox[0][2], 1.0,
        bbox[0][0], bbox[1][1], bbox[0][2], 1.0,
        bbox[0][0], bbox[0][1], bbox[1][2], 1.0,
        bbox[1][0], bbox[0][1], bbox[1][2], 1.0,
        bbox[1][0], bbox[1][1], bbox[1][2], 1.0,
        bbox[0][0], bbox[1][1], bbox[1][2], 1.0,
     };

    GLuint vbo_vertices;
    f->glGenBuffers(1, &vbo_vertices);
    f->glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
    f->glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    f->glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLushort elements[] = {
        0, 1, 2, 3,
        4, 5, 6, 7,
        0, 4, 1, 5, 2, 6, 3, 7
    };

    GLuint ibo_elements;
    f->glGenBuffers(1, &ibo_elements);
    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
    f->glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    f->glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
    int vertexLocation = program->attributeLocation("v_coord");
    program->enableAttributeArray(vertexLocation);
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 4, sizeof(float) * 4);

    f->glLineWidth(2.0);

    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
    f->glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, 0);
    f->glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, (GLvoid*)(4*sizeof(GLushort)));
    f->glDrawElements(GL_LINES, 8, GL_UNSIGNED_SHORT, (GLvoid*)(8*sizeof(GLushort)));
    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    program->disableAttributeArray(vertexLocation);
    f->glBindBuffer(GL_ARRAY_BUFFER, 0);

    f->glDeleteBuffers(1, &vbo_vertices);
    f->glDeleteBuffers(1, &ibo_elements);
}
