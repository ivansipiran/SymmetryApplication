! include(../common.pri){
    error("Could not add common config");
}

CONFIG += dll
CONFIG -= app_bundle
CONFIG -= qt

TARGET = Harris3D
TEMPLATE = lib

INCLUDEPATH += include/
INCLUDEPATH += ../Util/include
INCLUDEPATH += ../SimpleMesh/include

HEADERS += include/Harris3D/Distance.h \
           include/Harris3D/Face.h \
           include/Harris3D/HarrisDetector.h \
           include/Harris3D/Mesh.h \
           include/Harris3D/Vertex.h

SOURCES += src/Face.cpp \
           src/HarrisDetector.cpp \
           src/Mesh.cpp \
           src/Vertex.cpp

DEFINES += HARRIS_LIBRARY
LIBS += -LC:/Users/Admin/Documents/Sipiran/Implementations/tutorial_opengl/Qt/SymmetryApplication/bin/release -lSimpleMesh -lUtil
