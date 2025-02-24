#pragma once
#include "femcore/fem_export.h"

//-----------------------------------------------------------------------------
//! The FECore namespace encapsulates all classes that belong to the FECore library
namespace FECore
{
	// retrieve version numbers
	FEM_EXPORT void get_version(int& version, int& subversion);

	// retrieve version number string
	FEM_EXPORT const char* get_version_string();

	// initialize the module
	FEM_EXPORT void InitModule();

} // namespace FECore
