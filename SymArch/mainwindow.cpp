#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QtWidgets>
#include <QDir>
#include "SimpleMesh/io/offreader.h"
#include "SimpleMesh/io/offwriter.h"
#include "Harris3D/Mesh.h"
#include "Util/PropertySet.h"
#include "RotationalSymmetry/rotationalsymmetry.h"
#include "rotationaldialog.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    m_glArea = new GLWidget(this);
    setCentralWidget(m_glArea);
    createActions();
    createMenus();
    createToolbar();
    createStatusBar();
    createDockWindows();

    m_bOpenFile = false;
    setWindowTitle("SymArch");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::createDockWindows(){

    QDockWidget *dock = new QDockWidget(tr("Meshes"), this);
    dock->setFeatures(QDockWidget::NoDockWidgetFeatures);
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    m_ListMeshes = new QListWidget(dock);
    m_ListMeshes->setSelectionMode(QListWidget::ExtendedSelection);

    QWidget* dockContainer = new QWidget(this);
    QVBoxLayout* mainLayout = new QVBoxLayout;
    mainLayout->addWidget(m_ListMeshes);

    QPushButton* pushDeleteMesh = new QPushButton(QIcon("images/trash.png"), tr(""), this);
    QPushButton* pushMergeMeshes = new QPushButton(QIcon("images/merge.png"), tr(""), this);
    QPushButton* pushSegmentateMesh = new QPushButton(QIcon("images/segmentate.png"), tr(""), this);
    QPushButton* pushSaveMesh = new QPushButton(QIcon("images/save.png"), tr(""), this);

    connect(pushDeleteMesh, SIGNAL(pressed()), this, SLOT(deleteMesh()));
    connect(pushMergeMeshes, SIGNAL(pressed()), this, SLOT(mergeMeshes()));
    connect(pushSegmentateMesh, SIGNAL(pressed()), this, SLOT(segmentateMesh()));
    connect(pushSaveMesh, SIGNAL(pressed()), this, SLOT(saveMesh()));

    //QWidget* buttonContainer = new QWidget(this);
    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addWidget(pushDeleteMesh);
    buttonLayout->addWidget(pushMergeMeshes);
    buttonLayout->addWidget(pushSegmentateMesh);
    buttonLayout->addWidget(pushSaveMesh);
    buttonLayout->addStretch();
    //buttonContainer->setLayout(buttonLayout);
    mainLayout->addLayout(buttonLayout);

    dockContainer->setLayout(mainLayout);
    dock->setWidget(dockContainer);
    addDockWidget(Qt::RightDockWidgetArea, dock);

}

void MainWindow::createToolbar(){
    //m_toolbarMain = addToolBar(tr("Application"));
    //m_toolbarMain->addAction(m_actionOpen);
    ui->mainToolBar->addAction(m_actionOpen);
    m_toolbarMain = addToolBar(tr("Processing"));
    m_toolbarMain->addAction(m_actionRotSymmetry);
}

void MainWindow::createActions(){
    m_actionOpen = new QAction(QIcon("images/Open-icon.png"), tr("&Open"), this);
    m_actionOpen->setStatusTip(tr("Open a mesh"));
    connect(m_actionOpen, SIGNAL(triggered()), this, SLOT(open()));

    m_actionAbout = new QAction(tr("About"), this);
    m_actionAbout->setStatusTip(tr("About SymArch"));
    connect(m_actionAbout, SIGNAL(triggered()), this, SLOT(about()));

    m_actionExit = new QAction(tr("&Exit"), this);
    m_actionExit->setStatusTip(tr("Exit the application"));
    connect(m_actionExit, SIGNAL(triggered()), this, SLOT(close()));

    m_actionRotSymmetry = new QAction(QIcon("images/administrative-tools-icon.png"),tr("Rotational Symmetry"), this);
    m_actionRotSymmetry->setStatusTip("Repair rotationally symmetric objects");
    connect(m_actionRotSymmetry, SIGNAL(triggered()), this, SLOT(repairRotationalSymmetry()));
}

void MainWindow::createMenus(){
    m_menuFile = menuBar()->addMenu(tr("&File"));
    m_menuFile->addAction(m_actionOpen);
    m_menuFile->addSeparator();
    m_menuFile->addAction(m_actionExit);

    m_menuProcessing = menuBar()->addMenu(tr("&Processing"));
    m_menuProcessing->addAction(m_actionRotSymmetry);

    m_menuHelp = menuBar()->addMenu("&Help");
    m_menuHelp->addAction(m_actionAbout);
}

void MainWindow::createStatusBar(){
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::about(){
    QMessageBox::about(this, tr("About SymArch"),
                       tr("Symmetry Analysis and the application to archaeology"));
}

void MainWindow::open(){
    if(m_bOpenFile){ //Hay modelos abiertos, hay que gestionar la liberación de recursos
        m_glArea->clearContainer();
        m_bOpenFile = false;
        m_ListMeshes->clear();
    }

    QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Open mesh"), "./", tr("OFF Files (*.off)"));

    if(filenames.size()>0){
        QStringList list = filenames;
        QStringList::Iterator it = list.begin();
        while(it != list.end()){
            SimpleMesh::Mesh* new_mesh = new SimpleMesh::Mesh();
            SimpleMesh::IO::OFFReader reader(it->toUtf8().constData());
            reader.read_mesh(*new_mesh);
            new_mesh->compute_normals();

            SymArchData::BufferedObject3D* new_buffered_object = new SymArchData::BufferedObject3D();
            new_buffered_object->setGLWidget(m_glArea);
            new_buffered_object->createBuffersFromObject(new_mesh);
            int index = m_glArea->addObject(new_mesh, new_buffered_object);

            new_buffered_object->setColor(m_glArea->getColorPalette(index));

            //Extraer solo el nombre del archivo
            QStringList aux = it->split("/");
            qDebug() << QDir::separator();
            m_ListMeshes->addItem(aux.takeLast());
            ++it;

        }
        statusBar()->showMessage(tr("Object loaded"));
        //m_glArea->setMesh(m);
        m_bOpenFile = true;
        m_glArea->update();
    }
}

void MainWindow::saveMesh(){
    int index = m_ListMeshes->currentRow();

    if(index != -1){
        SimpleMesh::IO::OFFWriter writer("aaa.off");
        SimpleMesh::Mesh* m = m_glArea->getMesh(index);
        writer.write_mesh(*m);
    }
}

void MainWindow::deleteMesh(){
    int index = m_ListMeshes->currentRow();
    if(index!=-1){
       m_glArea->removeMesh(index);
       m_ListMeshes->takeItem(index);
       m_glArea->update();
    }
}

void MainWindow::mergeMeshes(){
    QModelIndexList indexes = m_ListMeshes->selectionModel()->selectedIndexes();
    std::vector<SimpleMesh::Mesh*> meshes;

    for(int i = 0; i < indexes.size(); i++){
        meshes.push_back(m_glArea->getMesh(indexes[i].row()));
    }

    std::cout << "List retrieved" << std::endl;

    if(meshes.size() > 1){ //Hacer la fusión sólo si hay más de una malla seleccionada
        unsigned int num_total_vertices = 0;
        unsigned int num_total_triangles = 0;

        for(unsigned int i = 0; i < meshes.size(); i++){
            num_total_vertices += meshes[i]->get_number_vertices();
            num_total_triangles += meshes[i]->get_number_triangles();
        }

        SimpleMesh::Mesh* bigMesh = new SimpleMesh::Mesh();
        bigMesh->set_number_vertices(num_total_vertices);
        bigMesh->set_number_triangles(num_total_triangles);

        std::cout << "Merged mesh created" << std::endl;

        unsigned int init_vertices = 0;
        unsigned int init_triangles = 0;

        for(unsigned int i = 0; i < meshes.size(); i++){
            SimpleMesh::Vertex* vertices = meshes[i]->get_vertices();
            SimpleMesh::Triangle* triangles = meshes[i]->get_triangles();

            for(unsigned int j = 0; j < meshes[i]->get_number_vertices(); j++)
                bigMesh->add_vertex(init_vertices + j, vertices[j].x(), vertices[j].y(), vertices[j].z());

            for(unsigned int  j = 0; j < meshes[i]->get_number_triangles(); j++)
                bigMesh->add_triangle(init_triangles + j, init_vertices + triangles[j].get_vertex_at(0), init_vertices + triangles[j].get_vertex_at(1), init_vertices + triangles[j].get_vertex_at(2));

            init_triangles += meshes[i]->get_number_triangles();
            init_vertices += meshes[i]->get_number_vertices();
        }

        std::cout << "Mesh merged" << std::endl;

        bigMesh->compute_normals();
        bigMesh->get_diagonal();

        std::cout << "Normals and diagonal" << std::endl;

        //Remover items de ListWidget y ObjectContainer
        for(int i = indexes.size()-1; i >= 0; i--){
            //std::cout << indexes[i].row() << std::endl;
            m_glArea->removeMesh(indexes[i].row());
            m_ListMeshes->takeItem(indexes[i].row());
        }

        std::cout << "Removed meshes and items" << std::endl;

        SymArchData::BufferedObject3D* new_buffered_object = new SymArchData::BufferedObject3D();
        new_buffered_object->setGLWidget(m_glArea);
        new_buffered_object->createBuffersFromObject(bigMesh);
        int ind = m_glArea->addObject(bigMesh, new_buffered_object);
        new_buffered_object->setColor(m_glArea->getColorPalette(ind));

        std::cout << "To ObjectContainer" << std::endl;

        m_ListMeshes->addItem("Merged");
        m_glArea->update();

    }



}

void MainWindow::repairRotationalSymmetry(){
    int index = m_ListMeshes->currentRow();

    if(index!=-1){
        //Leer parámetros desde diálogo
        Util::PropertySet prop;
        RotationalDialog* dialog = new RotationalDialog(this);
        dialog->setPropertySet(&prop);
        int result =  dialog->exec();

        if(result == QDialog::Accepted){
            SimpleMesh::Mesh* originalObject = m_glArea->getMesh(index);
            Harris3D::Mesh* copyObject = new Harris3D::Mesh();
            copyObject->deepCopy(originalObject);

            Engine::RotationalSymmetry engine(copyObject, &prop);
            SimpleMesh::Mesh* result = engine.computeCompletion();
            result->compute_normals();
            result->get_diagonal();

            SymArchData::BufferedObject3D* new_buffered_object = new SymArchData::BufferedObject3D();
            new_buffered_object->setGLWidget(m_glArea);
            new_buffered_object->createBuffersFromObject(result);
            int ind = m_glArea->addObject(result, new_buffered_object);
            new_buffered_object->setColor(m_glArea->getColorPalette(ind));

            m_ListMeshes->addItem("Result");
            m_glArea->update();
        }


        /*
        prop.addProperty(std::string("type-neighborhood"), std::string("rings"));
        prop.addProperty(std::string("parameter-neighborhood"), std::string("2"));
        prop.addProperty(std::string("K"), std::string("0.04"));
        prop.addProperty(std::string("ring-maxima-detection"), std::string("1"));
        prop.addProperty(std::string("interest-points-selection"), std::string("fraction"));
        prop.addProperty(std::string("parameter-selection"), std::string("0.01"));
        prop.addProperty(std::string("filtering-steps"), std::string("0"));
        prop.addProperty(std::string("alpha"), std::string("0.9"));
        prop.addProperty(std::string("factor-volumetric"), std::string("0.005"));
        prop.addProperty(std::string("num-triangles-decimation"), std::string("80000"));
        prop.addProperty(std::string("iter-non-rigid"), std::string("0"));
        prop.addProperty(std::string("factor-spacing"), std::string("2"));
        prop.addProperty(std::string("compatible-angle"), std::string("0.4"));
        prop.addProperty(std::string("enable-icp"), std::string("0"));
        prop.addProperty(std::string("rotation-factor"), std::string("2"));

        //Obtener la malla seleccionada
        SimpleMesh::Mesh* originalObject = m_glArea->getMesh(index);
        Harris3D::Mesh* copyObject = new Harris3D::Mesh();
        copyObject->deepCopy(originalObject);

        Engine::RotationalSymmetry engine(copyObject, &prop);
        SimpleMesh::Mesh* result = engine.computeCompletion();
        result->compute_normals();
        result->get_diagonal();

        SymArchData::BufferedObject3D* new_buffered_object = new SymArchData::BufferedObject3D();
        new_buffered_object->setGLWidget(m_glArea);
        new_buffered_object->createBuffersFromObject(result);
        int ind = m_glArea->addObject(result, new_buffered_object);
        new_buffered_object->setColor(m_glArea->getColorPalette(ind));

        m_ListMeshes->addItem("Result");
        m_glArea->update();*/

    }
}

void MainWindow :: segmentateMesh(){
    int index = m_ListMeshes->currentRow();

    if(index != -1){
        SimpleMesh::Mesh* object = m_glArea->getMesh(index);
        SimpleMesh::Vertex* vertices = object->get_vertices();
        SimpleMesh::Triangle* triangles = object->get_triangles();
        SimpleMesh::IO::OFFWriter writer("temp.off");
        writer.write_mesh(*object);

        Polyhedron p;
        std::ifstream input("temp.off");
        if(!input || !(input >> p) || p.empty()){
            std::cout << "Not a valid off file" << std::endl;
            return;
        }

        //Create a property map for sdf values
        typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
        Facet_double_map internal_sdf_map;
        boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

        //Compute sdf values
        CGAL::sdf_values(p, sdf_property_map);

        //Create a property map for segment-id
        typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
        Facet_int_map internal_segment_map;
        boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

        //Segmentate the mesh
        std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(p, sdf_property_map, segment_property_map);
        std::cout << "Number of segments:" << number_of_segments << endl;
        std::vector<SimpleMesh::Mesh*> meshes;

        for(unsigned int i = 0; i < number_of_segments; i++){
            //Identificar los vertices de este segmento
            std::vector<bool> mark;
            mark.resize(object->get_number_vertices());

            unsigned int num_faces = 0;
            for(Polyhedron::Facet_iterator facet_it = p.facets_begin(); facet_it!=p.facets_end(); ++facet_it){
                if(segment_property_map[facet_it] == i){
                    Polyhedron::Halfedge_around_facet_circulator circulator = facet_it->facet_begin();
                    int vertex0 = std::distance(p.vertices_begin(), circulator->vertex());  ++circulator;
                    int vertex1 = std::distance(p.vertices_begin(), circulator->vertex());  ++circulator;
                    int vertex2 = std::distance(p.vertices_begin(), circulator->vertex());
                    mark[vertex0] = true;   mark[vertex1] = true; mark[vertex2] = true;
                    num_faces++;
                }
            }

            unsigned int num_vertices = 0;

            std::map<unsigned int, unsigned int> mapping;
            for(unsigned int j = 0; j < mark.size(); j++){
                if(mark[j]){
                    mapping.insert(std::pair<unsigned int, unsigned int>(j, num_vertices));
                    num_vertices++;
                }
            }

            SimpleMesh::Mesh* auxMesh = new SimpleMesh::Mesh();
            auxMesh->set_number_vertices(num_vertices);
            auxMesh->set_number_triangles(num_faces);

            unsigned int cont = 0;
            for(unsigned int j = 0; j < mark.size(); j++){
                if(mark[j]){
                    auxMesh->add_vertex(cont++, vertices[j].x(), vertices[j].y(), vertices[j].z());
                }
            }

            cont = 0;
            for(Polyhedron::Facet_iterator facet_it = p.facets_begin(); facet_it!=p.facets_end(); ++facet_it){
                if(segment_property_map[facet_it] == i){
                    Polyhedron::Halfedge_around_facet_circulator circulator = facet_it->facet_begin();
                    int vertex0 = std::distance(p.vertices_begin(), circulator->vertex());  ++circulator;
                    int vertex1 = std::distance(p.vertices_begin(), circulator->vertex());  ++circulator;
                    int vertex2 = std::distance(p.vertices_begin(), circulator->vertex());
                    auxMesh->add_triangle(cont++, mapping[vertex0], mapping[vertex1], mapping[vertex2]);
                }
            }

            auxMesh->compute_normals();
            auxMesh->get_diagonal();
            if(auxMesh->get_number_vertices() > 100)
                meshes.push_back(auxMesh);
            else
                delete auxMesh;
        }

        for(int i = 0; i < meshes.size(); i++){
            SymArchData::BufferedObject3D* new_buffer = new SymArchData::BufferedObject3D();
            new_buffer->setGLWidget(m_glArea);
            new_buffer->createBuffersFromObject(meshes[i]);
            int ind = m_glArea->addObject(meshes[i], new_buffer);
            new_buffer->setColor(m_glArea->getColorPalette(ind));

            m_ListMeshes->addItem("Segment");
            m_glArea->update();
        }
    }
}
