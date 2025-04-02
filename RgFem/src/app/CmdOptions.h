/*********************************************************************
 * \file   CmdOptions.h
 * \brief  
 * 
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"
#include <string>

//! This structures stores the command line options that were input by the user
struct FEM_EXPORT CmdOptions
{
    constexpr static int MAXFILE = 512;

	int		ndebug;			//!< debug flag

	bool	bsplash;			//!< show splash screen or not
	bool	bsilent;			//!< run FEBio in silent mode (no output to screen)
	bool	binteractive;		//!< start FEBio interactively

	int		dumpLevel;		//!< requested restart level
	int		dumpStride;		//!< (cold) restart file stride

	char	szfile[MAXFILE];	//!< model input file name
	char	szlog[MAXFILE];	//!< log file name
	char	szplt[MAXFILE];	//!< plot file name
	char	szdmp[MAXFILE];	//!< dump file name
	char	szcnf[MAXFILE];	//!< configuration file
	char	sztask[MAXFILE];	//!< task name
	char	szctrl[MAXFILE];	//!< control file for tasks
	char	szimp[MAXFILE];		//!< import file

	CmdOptions()
	{
		defaults();
	}

	void defaults()
	{
		ndebug = 0;
		bsplash = true;
		bsilent = false;
		binteractive = false;
		dumpLevel = 0;
		dumpStride = 1;

		szfile[0] = 0;
		szlog[0] = 0;
		szplt[0] = 0;
		szdmp[0] = 0;
		szcnf[0] = 0;
		sztask[0] = 0;
		szctrl[0] = 0;
		szimp[0] = 0;
	}
};

// process string for command line options
FEM_EXPORT bool ProcessOptionsString(const std::string& s, CmdOptions& ops);
