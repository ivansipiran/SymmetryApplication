#ifndef UTIL_GLOBAL_H
#define UTIL_GLOBAL_H

//Helper macros to define library symbol visibility
#ifdef UTIL_EXPORT
#undef UTIL_EXPORT
#endif

#ifdef UTIL_IMPORT
#undef UTIL_IMPORT
#endif

#ifdef _MSC_VER
    #if defined(_DLL) && !defined(UTIL_STATICLIB) && !defined(UTIL_DLL)
        #define UTIL_DLL
    #endif
#endif

#ifdef __GNUC__
    #define UTIL_EXPORT __attribute__((visibility("default")))
    #define UTIL_IMPORT __attribute__((visibility("default")))
#endif

#ifdef _MSC_VER
    #ifdef UTIL_DLL
        #define UTIL_EXPORT __declspec(dllexport)
        #define UTIL_IMPORT __declspec(dllimport)
    #else
        #define UTIL_EXPORT
        #define UTIL_IMPORT
    #endif
#endif

#ifdef UTIL_API
#undef UTIL_API
#endif

#ifdef UTIL_LIBRARY
    #define UTIL_API UTIL_EXPORT
#else
    #define UTIL_API UTIL_IMPORT
#endif

#endif // UTIL_GLOBAL_H
