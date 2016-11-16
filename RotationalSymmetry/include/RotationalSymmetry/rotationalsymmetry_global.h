#ifndef ROTATIONALSYMMETRY_GLOBAL_H
#define ROTATIONALSYMMETRY_GLOBAL_H

//Helper macros to define library symbol visibility
#ifdef ROTATIONALSYMMETRY_EXPORT
#undef ROTATIONALSYMMETRY_EXPORT
#endif

#ifdef ROTATIONALSYMMETRY_IMPORT
#undef ROTATIONALSYMMETRY_IMPORT
#endif

#ifdef _MSC_VER
    #if defined(_DLL) && !defined(ROTATIONALSYMMETRY_STATICLIB) && !defined(ROTATIONALSYMMETRY_DLL)
        #define ROTATIONALSYMMETRY_DLL
    #endif
#endif

#ifdef __GNUC__
    #define ROTATIONALSYMMETRY_EXPORT __attribute__((visibility("default")))
    #define ROTATIONALSYMMETRY_IMPORT __attribute__((visibility("default")))
#endif

#ifdef _MSC_VER
    #ifdef ROTATIONALSYMMETRY_DLL
        #define ROTATIONALSYMMETRY_EXPORT __declspec(dllexport)
        #define ROTATIONALSYMMETRY_IMPORT __declspec(dllimport)
    #else
        #define ROTATIONALSYMMETRY_EXPORT
        #define ROTATIONALSYMMETRY_IMPORT
    #endif
#endif

#ifdef ROTATIONALSYMMETRY_API
#undef ROTATIONALSYMMETRY_API
#endif

#ifdef ROTATIONALSYMMETRY_LIBRARY
    #define ROTATIONALSYMMETRY_API ROTATIONALSYMMETRY_EXPORT
#else
    #define ROTATIONALSYMMETRY_API ROTATIONALSYMMETRY_IMPORT
#endif

#endif // SIMPLEMESH_GLOBAL_H

