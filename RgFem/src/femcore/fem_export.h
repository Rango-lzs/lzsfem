
#ifndef _FEM_EXPORT_H
#define _FEM_EXPORT_H

#define FEM_MODULE

#ifdef FEM_STATIC_DEFINE
#  define FEM_EXPORT
#  define FEM_NO_EXPORT
#else
#  ifndef FEM_EXPORT
#    ifdef FEM_MODULE
        /* We are building this library */
#      define FEM_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define FEM_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef FEM_NO_EXPORT
#    define FEM_NO_EXPORT 
#  endif
#endif

#ifndef FEM_DEPRECATED
#  define FEM_DEPRECATED __declspec(deprecated)
#endif

#ifndef FEM_DEPRECATED_EXPORT
#  define FEM_DEPRECATED_EXPORT FEM_EXPORT FEM_DEPRECATED
#endif

#ifndef FEM_DEPRECATED_NO_EXPORT
#  define FEM_DEPRECATED_NO_EXPORT FEM_NO_EXPORT FEM_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef FEM_NO_DEPRECATED
#    define FEM_NO_DEPRECATED
#  endif
#endif

#endif /* _FEM_EXPORT_H */
