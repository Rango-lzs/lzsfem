#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef fecore_EXPORTS
			#define FEM_EXPORT __declspec(dllexport)
		#else
			#define FEM_EXPORT __declspec(dllimport)
		#endif
	#else
		#define FEM_EXPORT
	#endif
#else
	#define FEM_EXPORT
#endif
