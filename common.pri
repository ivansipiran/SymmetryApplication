CONFIG *= warn_on ordered c++11  silent

unix|mingw|linux {
    CXXFLAGS *= -std=c11 -Wall -Wextra -Weffc++ -pedantic -Wno-parentheses
}

CONFIG (release, release|debug) {
    CXXFLAGS *= -O2
    unix|mingw {
        CXXFLAGS *= -mtune=generic -fomit-framepointer -ffast-math
    }
}

CONFIG (debug, debug|release) {
    DESTDIR = $${PWD}/bin/debug
    unix|mingw {
        CXXFLAGS *= -g -rdynamic -lmcheck
    }
}

CONFIG (release, release|debug) {
    DESTDIR = $${PWD}/bin/release
}

#LIBS += -L$${DESTDIR}

#-----------------CGAL------------------
unix {
    #message("Detected UNIX platform, using pkg-config to determine CGAL location")
    #CONFIG += link_pkgconfig
    #PKGCONFIG += cgal
    INCLUDEPATH += /usr/local/include
    LIBS += -L/usr/local/lib -lCGAL -lCGAL_Core
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/CGAL-release/include
    INCLUDEPATH += C:/Users/Admin/Documents/Software/CGAL-release/build/include
    LIBS += -LC:/Users/Admin/Documents/Software/CGAL-release/build/lib -lCGAL-vc120-mt-4.7 -lCGAL_Core-vc120-mt-4.7
}

#-----------------EIGEN------------------
unix {
    INCLUDEPATH += /usr/include/eigen3
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/eigen-eigen-07105f7124f9
}

#-----------------HALF------------------
unix {
    INCLUDEPATH += /usr/include
    LIBS += -L/usr/lib/x86_64-linux-gnu -lHalf -lhdf5
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/ilmbase-2.2.0/ilmbase-2.2.0/Half
    LIBS += -LC:/Users/Admin/Documents/Software/ilmbase-2.2.0/ilmbase-2.2.0/build/Half/Release -lHalf
}

#-----------------TBB------------------
unix {
    INCLUDEPATH += /home/sipiran/Downloads/tbb44_20160128oss/include
    LIBS += -L/home/sipiran/Downloads/tbb44_20160128oss/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.16.0_release -ltbb
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/tbb44_20160128oss_src/tbb44_20160128oss/include
    LIBS += -LC:/Users/Admin/Documents/Software/tbb44_20160128oss_src/tbb44_20160128oss/build/vs2010/x64/Release -ltbb
}

#-----------------LOG4CPLUS------------------
unix {
    INCLUDEPATH += /usr/include/log4cplus
    LIBS += -llog4cplus
}else{

}

#-----------------OPENVDB------------------
unix {
    INCLUDEPATH += /home/sipiran/Downloads
    LIBS += -L/home/sipiran/Downloads/openvdb -lopenvdb
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/openvdb_3_1_0_library
    #LIBS += C:/Users/Admin/Documents/Software/openvdb_3_1_0_library/openvdb/OpenVDB.lib
}

#-----------------BOOST------------------
unix {
    INCLUDEPATH += /usr/include
    LIBS += -L/usr/lib/x86_64-linux-gnu -lboost_system -lboost_iostreams -llog4cplus -lboost_thread

}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/boost-release/include/boost-1_60
    LIBS += -LC:/Users/Admin/Documents/Software/boost-release/lib -lboost_thread-vc120-mt-1_60
}

#----------------PCL---------------------
unix{
    message("Detected UNIX platform, using pkg-config to determine PCL location")
    INCLUDEPATH += /usr/local/include/pcl-1.7
    LIBS += -lpcl_common -lpcl_kdtree -lpcl_sample_consensus -lpcl_search
    #PKGCONFIG += pcl_common-1.7 pcl_kdtree-1.7 pcl_sample_consensus-1.7 pcl_search-1.7
}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/pcl-pcl-1.7.2/include
    LIBS += -LC:/Users/Admin/Documents/Software/pcl-pcl-1.7.2/build/lib -lpcl_common_release -lpcl_kdtree_release -lpcl_sample_consensus_release -lpcl_search_release
}

#-----------------FLANN-------------------
unix{

}else{
    INCLUDEPATH += C:/Users/Admin/Documents/Software/flann-1.8.4-src/src/cpp
    LIBS += -LC:/Users/Admin/Documents/Software/flann-1.8.4-src/build/lib/Release
}
