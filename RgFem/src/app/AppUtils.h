#pragma once
#include <cstddef>
#include "femcore/fem_export.h"
#include "FEAppConfig.h"
#include "cmdoptions.h"
#include <vector>
#include <ostream>

class CompactMatrix;
class LogStream;
class FeModel;
class FEMaterial;

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace Rango 
{
	FEM_EXPORT FEModel* CreateFEModel();
	// Initialize all the FEBio modules
	FEM_EXPORT void InitLibrary();

	// read the configuration file
	FEM_EXPORT bool Configure(const char* szfile, FEAppConfig& config);

	// load a plugin
	FEM_EXPORT bool ImportPlugin(const char* szfile);

	// load all the plugins in a folder
	FEM_EXPORT void ImportPluginFolder(const char* szfolder);

	// get the name of the plugin from its allocator Id
	FEM_EXPORT const char* GetPluginName(int allocId);

	// call this to clean up all FEBio data
	FEM_EXPORT void FinishLibrary();

	// helper function for retrieving the executable's path
	FEM_EXPORT int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FEM_EXPORT int Hello(LogStream& log);

	// set the number of OMP threads
	FEM_EXPORT void SetOMPThreads(int n);

	// run an FEBioModel
	FEM_EXPORT bool SolveModel(FeModel& fem, const char* sztask = nullptr, const char* szctrl = nullptr);

	// run an FEBioModel
	FEM_EXPORT int RunModel(FeModel& fem, CmdOptions* ops);

	// write a matrix to file
	FEM_EXPORT bool write_hb(CompactMatrix& K, const char* szfile, int mode = 0);

	// print matrix sparsity pattern to svn file
	FEM_EXPORT void print_svg(CompactMatrix* m, std::ostream &out, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

	// write a vector to file
	FEM_EXPORT bool write_vector(const std::vector<double>& a, const char* szfile, int mode = 0);

	// run a material test
	FEM_EXPORT bool RunMaterialTest(FEMaterial* mat, double simtime, int steps, double strain, const char* sztest, std::vector<std::pair<double, double> >& out);
}
