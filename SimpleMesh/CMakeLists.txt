cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project ( SimpleMesh )

set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/io/objreader.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/io/offreader.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/io/objwriter.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/io/offwriter.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/mesh.cpp )

add_definitions( -DSIMPLEMESH_LIBRARY )
add_definitions( -DDEBUG )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/.. )
include_directories( C:/Users/Admin/Documents/Software/boost-release/include/boost-1_60 )
add_library( SimpleMesh SHARED ${SOURCES})
