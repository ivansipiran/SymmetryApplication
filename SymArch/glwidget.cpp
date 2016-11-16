#include "glwidget.h"
#include <QMouseEvent>
#include <QOpenGLShaderProgram>
#include <QCoreApplication>
#include <QtMath>
#include <math.h>
#include <iostream>

#include "SimpleMesh/io/offreader.h"
#include "SimpleMesh/io/offwriter.h"
#include "analysisbar.h"


using namespace std;
GLWidget::GLWidget(QWidget* parent)
    :QOpenGLWidget(parent),
      ArcBall(400.0f, 400.0f)

{
    m_core = QCoreApplication::arguments().contains(QStringLiteral("--coreprofile"));
    // --transparent causes the clear color to be transparent. Therefore, on systems that
    // support it, the widget will become transparent apart from the logo.
    m_transparent = QCoreApplication::arguments().contains(QStringLiteral("--transparent"));
    if (m_transparent)
        setAttribute(Qt::WA_TranslucentBackground);

    Transform.setToIdentity();
    LastRot.setToIdentity();
    ThisRot.setToIdentity();
    isClicked = false;
    delta_scale = 0.0;

    m_programs[MAIN_OBJECT_PROGRAM] = NULL;
    m_programs[BOUNDING_BOX_PROGRAM] = NULL;
    m_programs[COLORBAR_PROGRAM] = NULL;


    //Inicializar la paleta de colores
    palette[0].setX(0.86);  palette[0].setY(0.86);  palette[0].setZ(0.86); //Gris claro, color por defecto
    palette[1].setX(0.84);  palette[1].setY(1.00);  palette[1].setZ(0.24); //Verde limon
    palette[2].setX(1.00);  palette[2].setY(0.80);  palette[2].setZ(0.40); //Maiz
    palette[3].setX(0.80);  palette[3].setY(0.90);  palette[3].setZ(1.00); //Celeste
    palette[4].setX(0.79);  palette[4].setY(0.58);  palette[4].setZ(0.89); //Morado
    palette[5].setX(0.00);  palette[5].setY(0.60);  palette[5].setZ(1.00); //Azul
    palette[6].setX(1.00);  palette[6].setY(1.00);  palette[6].setZ(0.24); //Amarillo
    palette[7].setX(0.60);  palette[7].setY(0.60);  palette[7].setZ(0.60); //Gris oscuro


    //palette[8].setX(0.86);  palette[8].setY(0.86);  palette[8].setZ(0.86); //
    //palette[9].setX(0.86);  palette[9].setY(0.86);  palette[9].setZ(0.86); //
    //palette[10].setX(0.86);  palette[10].setY(0.86);  palette[10].setZ(0.86); //
    //palette[11].setX(0.86);  palette[11].setY(0.86);  palette[11].setZ(0.86); //
    //palette[12].setX(0.86);  palette[12].setY(0.86);  palette[12].setZ(0.86); //
    //palette[13].setX(0.86);  palette[13].setY(0.86);  palette[13].setZ(0.86); //
    //palette[14].setX(0.86);  palette[14].setY(0.86);  palette[14].setZ(0.86); //
    //palette[15].setX(0.86);  palette[15].setY(0.86);  palette[15].setZ(0.86); //
}

GLWidget::~GLWidget(){
    cleanup();
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}

static void qNormalizeAngle(int &angle){
    while(angle < 0)
        angle += 360*16;
    while(angle > 360*16)
        angle -=360*16;
}

void GLWidget::cleanup(){
    makeCurrent();
    for(int i = 0; i < NUM_PROGRAMS; i++)
        if(m_programs[i]!=NULL)
            delete m_programs[i];
    doneCurrent();
}

int GLWidget::addObject(SimpleMesh::Mesh* m, SymArchData::BufferedObject3D* bo){
     int size = m_container.getNumberObjects();
     m_container.addObject(m, bo);
     return size;
}

void GLWidget::initializeGL(){
    bool ok = true;
    connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);

    initializeOpenGLFunctions();
    glClearColor(0, 0, 0, m_transparent ? 0 : 1);

    //Set the shaders for main object
    m_programs[MAIN_OBJECT_PROGRAM] = new QOpenGLShaderProgram;
    ok = m_programs[MAIN_OBJECT_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Vertex, "phong-shading.v.glsl");
    if(!ok){
        cout << "Problem compiling Vertex Shader" << endl;
        cout << m_programs[MAIN_OBJECT_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }

    ok = m_programs[MAIN_OBJECT_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Fragment, "phong-shading.f.glsl");
    if(!ok){
        cout << "Problem compiling Fragment Shader"<< endl;
        cout << m_programs[MAIN_OBJECT_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }


    m_programs[MAIN_OBJECT_PROGRAM]->link();

    //Set the shaders for bounding box
    m_programs[BOUNDING_BOX_PROGRAM] = new QOpenGLShaderProgram;
    ok = m_programs[BOUNDING_BOX_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Vertex, "bounding-shading.v.glsl");
    if(!ok){
        cout << "Problem compiling Vertex Shader" << endl;
        cout << m_programs[BOUNDING_BOX_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }

    ok = m_programs[BOUNDING_BOX_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Fragment, "bounding-shading.f.glsl");
    if(!ok){
        cout << "Problem compiling Fragment Shader" << endl;
        cout << m_programs[BOUNDING_BOX_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }

    m_programs[BOUNDING_BOX_PROGRAM]->link();

    //Set the shaders for per-vertex color objects
   /* m_programs[COLORBAR_PROGRAM] = new QOpenGLShaderProgram;
    ok = m_programs[COLORBAR_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Vertex, "texture.v.glsl");
    if(!ok){
        cout << "Problem compiling Vertex Shader" << endl;
        cout << m_programs[COLORBAR_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }

    ok = m_programs[COLORBAR_PROGRAM]->addShaderFromSourceFile(QOpenGLShader::Fragment, "texture.f.glsl");
    if(!ok){
        cout << "Problem compiling Fragment Shader" << endl;
        cout << m_programs[COLORBAR_PROGRAM]->log().toStdString() << endl;
        qApp->exit(-1);
    }

    m_programs[COLORBAR_PROGRAM]->link();*/

    //m_programs[COLOR_PER_VERTEX_PROGRAM]->bind();
        //attribute_v_coord = m_programs[COLOR_PER_VERTEX_PROGRAM]->attributeLocation("v_coord");
        //attribute_v_normal = m_programs[COLOR_PER_VERTEX_PROGRAM]->attributeLocation("v_normal");
    //    uniform_m = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("m");
    //    uniform_v = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("v");
    //    uniform_p = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("p");
    //    uniform_m_3x3_inv_transp = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("m_3x3_inv_transp");
    //    uniform_v_inv = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("v_inv");
    //m_programs[COLOR_PER_VERTEX_PROGRAM]->release();

    m_programs[MAIN_OBJECT_PROGRAM]->bind();
        uniform_m = m_programs[MAIN_OBJECT_PROGRAM]->uniformLocation("m");
        uniform_v = m_programs[MAIN_OBJECT_PROGRAM]->uniformLocation("v");
        uniform_p = m_programs[MAIN_OBJECT_PROGRAM]->uniformLocation("p");
        uniform_m_3x3_inv_transp = m_programs[MAIN_OBJECT_PROGRAM]->uniformLocation("m_3x3_inv_transp");
        uniform_v_inv = m_programs[MAIN_OBJECT_PROGRAM]->uniformLocation("v_inv");
    m_programs[MAIN_OBJECT_PROGRAM]->release();

    //SimpleMesh::Mesh* new_mesh = new SimpleMesh::Mesh();
    //SimpleMesh::IO::OFFReader reader("814_Ramesses_1.5Mtriangles_clean.off");
    //reader.read_mesh(*new_mesh);
    //new_mesh->compute_normals();

    //AnalysisBar* new_mesh = new AnalysisBar(64, 64, 0.05);
    //new_mesh->set_flag_color(true);
    //new_mesh->proceduralInit();

    //SimpleMesh::IO::OFFWriter writer("a.off");
    //writer.write_mesh(*new_mesh);

    //SymArchData::BufferedObject3D* new_buffered_object = new SymArchData::BufferedObject3D();
    //new_buffered_object->createBuffersFromObject(new_mesh);
    //m_container.addObject(new_mesh, new_buffered_object);

    m_cam.setToIdentity();
    m_cam.translate(0, 0, -1);
    m_cam.lookAt(QVector3D(0.0, 0.0, 4.0), QVector3D(0.0, 0.0, 0.0), QVector3D(0.0, 1.0, 0.0));



}

void GLWidget::paintGL(){

    scale = 1.5*m_container.getScale();
    center = m_container.getCenter();

    QMatrix4x4 translation;
    translation.setToIdentity();
    translation.translate(-center);

    QMatrix4x4 scal;
    scal.setToIdentity();
    scal.scale(scale);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    m_programs[MAIN_OBJECT_PROGRAM]->bind();
        m_programs[MAIN_OBJECT_PROGRAM]->setUniformValue(uniform_v, m_cam);
        m_programs[MAIN_OBJECT_PROGRAM]->setUniformValue(uniform_p, m_proj);
        QMatrix4x4 v_inv = m_cam.inverted();
        m_programs[MAIN_OBJECT_PROGRAM]->setUniformValue(uniform_v_inv, v_inv);
        m_programs[MAIN_OBJECT_PROGRAM]->setUniformValue(uniform_m, Transform*scal*translation);
        QMatrix3x3 m_3x3_inv_transp = Transform.normalMatrix();
        m_programs[MAIN_OBJECT_PROGRAM]->setUniformValue(uniform_m_3x3_inv_transp, m_3x3_inv_transp);

        m_container.drawObjects(m_programs[MAIN_OBJECT_PROGRAM]);
     m_programs[MAIN_OBJECT_PROGRAM]->release();

    //m_programs[COLOR_PER_VERTEX_PROGRAM]->bind();
    //    uniform_m = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("m");
    //    uniform_v = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("v");
    //    uniform_p = m_programs[COLOR_PER_VERTEX_PROGRAM]->uniformLocation("p");

    //    m_programs[COLOR_PER_VERTEX_PROGRAM]->setUniformValue(uniform_v, m_cam);
    //    m_programs[COLOR_PER_VERTEX_PROGRAM]->setUniformValue(uniform_p, m_proj);
    //    QMatrix4x4 v_inv = m_cam.inverted();
    //    m_programs[COLOR_PER_VERTEX_PROGRAM]->setUniformValue(uniform_v_inv, v_inv);
    //    m_programs[COLOR_PER_VERTEX_PROGRAM]->setUniformValue(uniform_m, Transform*scal*translation);
    //    QMatrix3x3 m_3x3_inv_transp = Transform.normalMatrix();
    //    m_programs[COLOR_PER_VERTEX_PROGRAM]->setUniformValue(uniform_m_3x3_inv_transp, m_3x3_inv_transp);
    //    m_container.drawObjects(m_programs[COLOR_PER_VERTEX_PROGRAM]);
    // m_programs[COLOR_PER_VERTEX_PROGRAM]->release();


     //m_programs[BOUNDING_BOX_PROGRAM]->bind();
     //   m_programs[BOUNDING_BOX_PROGRAM]->setUniformValue("v", m_cam);
     //   m_programs[BOUNDING_BOX_PROGRAM]->setUniformValue("m", Transform*scal*translation);
     //   m_programs[BOUNDING_BOX_PROGRAM]->setUniformValue("p", m_proj);
     //   m_container.drawBoundingBox(m_programs[BOUNDING_BOX_PROGRAM]);
     //m_programs[BOUNDING_BOX_PROGRAM]->release();


    glFlush();


}

void GLWidget::resizeGL(int width, int height){
    screen_width = width;
    screen_height = height;
    m_proj.setToIdentity();
    m_proj.perspective(45.0f, GLfloat(width) / height, 0.01f, 100.0f);
    ArcBall.setBounds((GLfloat)width, (GLfloat)height);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    isClicked = true;
    LastRot = ThisRot;
    Point2fT auxPoint;
    auxPoint.s.X = (GLfloat)event->pos().x();
    auxPoint.s.Y = (GLfloat)event->pos().y();
    ArcBall.click(&auxPoint);

}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if(isClicked){
        Quat4fT ThisQuat;
        Point2fT auxPoint;
        auxPoint.s.X = (GLfloat)event->pos().x();
        auxPoint.s.Y = (GLfloat)event->pos().y();
        ArcBall.drag(&auxPoint, &ThisQuat);
        Matrix3fT mThisRot;
        memcpy(mThisRot.M, ThisRot.data(), sizeof(float)*9);

        Matrix3fT mLastRot;
        memcpy(mLastRot.M, LastRot.data(), sizeof(float)*9);
        Matrix4fT mTransform;
        memcpy(mTransform.M, Transform.data(), sizeof(float)*16);

        Matrix3fSetRotationFromQuat4f(&mThisRot, &ThisQuat);
        Matrix3fMulMatrix3f(&mThisRot, &mLastRot);
        Matrix4fSetRotationFromMatrix3f(&mTransform, &mThisRot);

        memcpy(ThisRot.data(), mThisRot.M, sizeof(float)*9);
        memcpy(LastRot.data(), mLastRot.M, sizeof(float)*9);
        memcpy(Transform.data(), mTransform.M, sizeof(float)*16);
        update();
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    isClicked = false;
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    //float delta = (float)event->delta();
    //delta_scale = delta/12000.0;
    //update();
}
