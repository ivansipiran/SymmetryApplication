#ifndef HARRIS_GLOBAL_H
#define HARRIS_GLOBAL_H

//Helper macros to define library symbol visibility
#ifdef HARRIS_EXPORT
#undef HARRIS_EXPORT
#endif

#ifdef HARRIS_IMPORT
#undef HARRIS_IMPORT
#endif

#ifdef _MSC_VER
#if defined(_DLL) && !defined(HARRIS_STATICLIB) && !defined(HARRIS_DLL)
#define HARRIS_DLL
#endif
#endif

#ifdef __GNUC__
#define HARRIS_EXPORT __attribute__((visibility("default")))
#define HARRIS_IMPORT __attribute__((visibility("default")))
#endif

#ifdef _MSC_VER
#ifdef HARRIS_DLL
#define HARRIS_EXPORT __declspec(dllexport)
#define HARRIS_IMPORT __declspec(dllimport)
#else
#define HARRIS_EXPORT
#define HARRIS_IMPORT
#endif
#endif

#ifdef HARRIS_API
#undef HARRIS_API
#endif

#ifdef HARRIS_LIBRARY
#define HARRIS_API HARRIS_EXPORT
#else
#define HARRIS_API HARRIS_IMPORT
#endif

#endif // HARRIS_GLOBAL_H
