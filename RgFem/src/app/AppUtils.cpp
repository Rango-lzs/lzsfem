#include "app/AppUtils.h"
#include "femcore/FERgModel.h"

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace Rango 
{
	FEModel* CreateFEModel()
	{
        return new FERgModel();
	}

	// Initialize all the FEBio modules
	void InitLibrary()
    {
    }

	// read the configuration file
	bool Configure(const char* szfile, FEAppConfig& config)
    {
        return true;
    }

	// load a plugin
	bool ImportPlugin(const char* szfile)
    {
        return true;
    }

	// load all the plugins in a folder
    void ImportPluginFolder(const char* szfolder)
    {
    }

	// get the name of the plugin from its allocator Id
	const char* GetPluginName(int allocId)
    {
        return nullptr;
    }

	// call this to clean up all FEBio data
	void FinishLibrary()
    {
    }

	// helper function for retrieving the executable's path
	int get_app_path(char *pname, size_t pathsize)
    {
        return 0;
	}

	// print hello message
	int Hello(LogStream& log)
    {
        return 0;
    }

	// set the number of OMP threads
	void SetOMPThreads(int n)
    {
        
    }

	// run an FEBioModel
	bool SolveModel(FEModel& fem, const char* sztask, const char* szctrl)
    {
        return true;
    }

	// run an FEBioModel
	int RunModel(FEModel& fem, CmdOptions* ops)
    {
        return 0;
    }

	// write a Matrix to file
	bool write_hb(CompactMatrix& K, const char* szfile, int mode)
    {
        return true;
    }

	// print Matrix sparsity pattern to svn file
	void print_svg(CompactMatrix* m, std::ostream &out, int i0, int j0, int i1, int j1)
    {
    }

	// write a vector to file
	bool write_vector(const std::vector<double>& a, const char* szfile, int mode)
    {
        return true;
    }

	// run a material test
	bool RunMaterialTest(FEMaterial* mat, double simtime, int steps, double strain, const char* sztest, std::vector<std::pair<double, double> >& out)
    {
        return true;
    }
}
