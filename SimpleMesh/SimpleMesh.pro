! include(../common.pri){
    error("Could not add common config");
}

CONFIG += dll
CONFIG -= app_bundle
CONFIG -= qt

TARGET = SimpleMesh
TEMPLATE = lib

INCLUDEPATH += include/

HEADERS += include/SimpleMesh/simplemesh_global.h \
           include/SimpleMesh/vertex.h \
           include/SimpleMesh/triangle.h \
           include/SimpleMesh/mesh.h \
           include/SimpleMesh/io/offreader.h \
           include/SimpleMesh/io/objreader.h \
           include/SimpleMesh/io/offwriter.h \
           include/SimpleMesh/io/objwriter.h

SOURCES +=  src/SimpleMesh/mesh.cpp \
            src/SimpleMesh/io/offreader.cpp \
            src/SimpleMesh/io/objreader.cpp \
            src/SimpleMesh/io/offwriter.cpp \
            src/SimpleMesh/io/objwriter.cpp

DEFINES += SIMPLEMESH_LIBRARY
