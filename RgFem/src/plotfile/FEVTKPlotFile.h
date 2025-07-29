#pragma once
#include "PlotFile.h"
#include <fstream>

//-----------------------------------------------------------------------------
//! This class implements the facilities to export FE data in the VTK
//! plot file format.
//!
class FEVTKPlotFile : public PlotFile
{
public:
    FEVTKPlotFile(FEModel* fem)
        : PlotFile(fem)
    {
    }

	~FEVTKPlotFile(void)
    {
    }

	//! Open the plot database
	bool Open(const char* szfile) override;

	//! Close the plot database
	void Close() override;

	//! Open for appending
	bool Append(const char* szfile) override;

	//! Write current FE state to plot database
	bool Write(float ftime, int flag = 0)  override;

	//! see if the plot file is valid
	bool IsValid() const override;

private:
    std::ofstream m_out;
};
