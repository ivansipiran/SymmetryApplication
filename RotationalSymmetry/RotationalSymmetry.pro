! include(../common.pri){
    error("Could not add common config");
}

QT       -= core gui

TARGET = RotationalSymmetry
TEMPLATE = lib

DEFINES += ROTATIONALSYMMETRY_LIBRARY

SOURCES += src/rotationalsymmetry.cpp

HEADERS += include/RotationalSymmetry/rotationalsymmetry.h\
        include/RotationalSymmetry/rotationalsymmetry_global.h

INCLUDEPATH += ./include
INCLUDEPATH += ../Util/include
INCLUDEPATH += ../Harris3D/include
INCLUDEPATH += ../SimpleMesh/include
INCLUDEPATH += ../sparseicp
INCLUDEPATH += ../sparseicp/nanoflann-1.1.7/include

LIBS += -LC:/Users/Admin/Documents/Sipiran/Implementations/tutorial_opengl/Qt/SymmetryApplication/bin/release -lSimpleMesh -lUtil -lHarris3D
LIBS += -LC:/Users/Admin/Documents/Software/openvdb_3_1_0_library/openvdb/build -lOpenVDB
