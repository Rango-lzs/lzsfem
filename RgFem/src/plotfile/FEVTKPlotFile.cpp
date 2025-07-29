#include "FEVTKPlotFile.h"
#include <iostream>
#include "VtuWriter.h"
#include "femcore/FEModel.h"

#include <vector>

//-----------------------------------------------------------------------------
bool FEVTKPlotFile::IsValid() const
{
	//return m_ar.IsValid();
    return false;
}

//-----------------------------------------------------------------------------
void FEVTKPlotFile::Close()
{
    m_out.close();
}


//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Open(const char *szfile)
{
	FEModel* fem = GetFEModel();
    if (!fem)
    {
        std::cerr << "Error: No FEModel available." << std::endl;
        return false;
    }

    m_out.open(szfile);
    if (!m_out.is_open())
    {
        std::cerr << "Error: Could not open file " << szfile << std::endl;
        return false;
    }   
	return true;
}

//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Write(float ftime, int flag)
{
	FEModel& fem = *GetFEModel();
    VTUWriter vtuOut(&fem, m_out);
    vtuOut.write(std::vector<NodeData>(), std::vector<ElementData>());
	return true;
}

//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Append(const char *szfile)
{
	return false;
}
