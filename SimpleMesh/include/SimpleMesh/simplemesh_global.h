#ifndef SIMPLEMESH_GLOBAL_H
#define SIMPLEMESH_GLOBAL_H

//Helper macros to define library symbol visibility
#ifdef SIMPLEMESH_EXPORT
#undef SIMPLEMESH_EXPORT
#endif

#ifdef SIMPLEMESH_IMPORT
#undef SIMPLEMESH_IMPORT
#endif

#ifdef _MSC_VER
    #if defined(_DLL) && !defined(SIMPLEMESH_STATICLIB) && !defined(SIMPLEMESH_DLL)
        #define SIMPLEMESH_DLL
    #endif
#endif

#ifdef __GNUC__
    #define SIMPLEMESH_EXPORT __attribute__((visibility("default")))
    #define SIMPLEMESH_IMPORT __attribute__((visibility("default")))
#endif

#ifdef _MSC_VER
    #ifdef SIMPLEMESH_DLL
        #define SIMPLEMESH_EXPORT __declspec(dllexport)
        #define SIMPLEMESH_IMPORT __declspec(dllimport)
    #else
        #define SIMPLEMESH_EXPORT
        #define SIMPLEMESH_IMPORT
    #endif
#endif

#ifdef SIMPLEMESH_API
#undef SIMPLEMESH_API
#endif

#ifdef SIMPLEMESH_LIBRARY
    #define SIMPLEMESH_API SIMPLEMESH_EXPORT
#else
    #define SIMPLEMESH_API SIMPLEMESH_IMPORT
#endif

#endif // SIMPLEMESH_GLOBAL_H
