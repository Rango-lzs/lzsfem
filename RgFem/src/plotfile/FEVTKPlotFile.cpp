#include "FEVTKPlotFile.h"


//-----------------------------------------------------------------------------
bool FEVTKPlotFile::IsValid() const
{
	//return m_ar.IsValid();
    return false;
}

//-----------------------------------------------------------------------------
void FEVTKPlotFile::Close()
{
	//m_ar.Close();
}


//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Open(const char *szfile)
{
	FEModel* fem = GetFEModel();

	// open the archive
	//m_ar.Create(szfile);

	return true;
}

//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Write(float ftime, int flag)
{
	FEModel& fem = *GetFEModel();
	return true;
}

//-----------------------------------------------------------------------------
bool FEVTKPlotFile::Append(const char *szfile)
{
	return false;
}
