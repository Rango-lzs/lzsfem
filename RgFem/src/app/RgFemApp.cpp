/*********************************************************************
 * \file   RgFemApp.cpp
 * \brief  
 * 
 * \author Leizs
 * \date   December 2024
 *********************************************************************/

#include "RgFemApp.h"
#include "app/CmdOptions.h"

#include "CLI/CLI.hpp"

RgFemApp::RgFemApp()
{
    mp_fem = nullptr;
}

RgFemApp* RgFemApp::Instance()
{
	static RgFemApp theApp;
	return &theApp;
}

bool RgFemApp::Init(int argc, char* argv[])
{

	// 模块的初始化,类的注册,类的信息可以在模块加载时自动注册，也可以手动去初始化~
	// 这部分等一些基础类写好后，再进行注册模块的注册
	//FeModuleInit::InitLibrary();

	CLI::App femApp("Rango FEM Application");
    std::string inputFile;
    femApp.add_option("-f,-file", inputFile, "the fem input file");
	femApp.parse(argc,argv);

	ParseCmdLine(argc, argv);

	// copy some flags to configuration
	m_config.SetOutputLevel(m_cmd_opts->bsilent ? 0 : 1);

	// read the configration file if specified
	if (m_cmd_opts.szcnf[0])
		if (ReadConfigure(m_cmd_opts.szcnf, m_config) == false)
		{
			fprintf(stderr, "FATAL ERROR: An error occurred reading the configuration file.\n");
			return false;
		}

	// read command line plugin if specified, 可以使用开源库来代替,可以先不实现这部分
	if (m_cmd_opts.szimp[0] != 0)
	{
		ImportPlugin(m_cmd_opts.szimp);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool RgFemApp::ReadConfigure(const char* szconfig)
{
	return Configure(szconfig, m_config);
}

//-----------------------------------------------------------------------------
int RgFemApp::Run()
{
	return RunModel();
}

//-----------------------------------------------------------------------------
void RgFemApp::Finish()
{
	febio::FinishLibrary();

	Console::GetHandle()->CleanUp();
}

//-----------------------------------------------------------------------------
// get the current model
FEBioModel* RgFemApp::GetCurrentModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
// set the currently active model
void RgFemApp::SetCurrentModel(FEBioModel* fem)
{
	m_fem = fem;
}

//-----------------------------------------------------------------------------
// Run an input file. 
int RgFemApp::RunModel()
{
	// create the FEModel object
	FEBioModel model;
    SetCurrentModel(&model);

	// read the input file if specified
	if (m_cont.inPutFile[0])
	{
		// read the input file
        if (!model.Input(m_cmd_opts.szfile))
		{
			return 1;
		}

	    // apply configuration overrides
        ApplyConfig(model);
	}

	// solve the model with the task and control file
    bool ret = febio::SolveModel(model, m_cmd_opts.sztask, m_cmd_opts.szctrl);
	
	// reset the current model pointer
	SetCurrentModel(nullptr);
	return  ret? 1 : 0;
}

//-----------------------------------------------------------------------------
// apply configuration changes to model
void RgFemApp::ApplyConfig(FEBioModel& fem)
{
	if (mp_config.m_printParams != -1)
	{
		fem.SetPrintParametersFlag(mp_config.m_printParams != 0);
	}
	fem.ShowWarningsAndErrors(mp_config.m_bshowErrors);
}

//-----------------------------------------------------------------------------
//!  Parses the command line and returns a CMDOPTIONS structure
//
bool RgFemApp::ParseCmdLine(int nargs, char* argv[])
{
	CmdOptions& ops = *mp_cmd_opts;

	// set default options
	ops.ndebug = 0;
	ops.bsplash = true;
	ops.bsilent = false;
	ops.binteractive = true;

	// these flags indicate whether the corresponding file name
	// was defined on the command line. Otherwise, a default name will be generated.
	bool blog = false;
	bool bplt = false;
	bool bdmp = false;
	bool brun = true;

	// initialize file names
	ops.szfile[0] = 0;
	ops.szplt[0] = 0;
	ops.szlog[0] = 0;
	ops.szdmp[0] = 0;
	ops.sztask[0] = 0;
	ops.szctrl[0] = 0;
	ops.szimp[0] = 0;

	// set initial configuration file name
	if (ops.szcnf[0] == 0)
	{
		char szpath[1024] = { 0 };
		febio::get_app_path(szpath, 1023);
		sprintf(ops.szcnf, "%sfebio.xml", szpath);
	}

	// loop over the arguments
	char* sz;
	for (int i=1; i<nargs; ++i)
	{
		sz = argv[i];

		if (strcmp(sz,"-r") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-r is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "restart");
			strcpy(ops.szctrl, argv[++i]);
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-d") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-d is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "diagnose");
			strcpy(ops.szctrl, argv[++i]);
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-p") == 0)
		{
			bplt = true;
			strcpy(ops.szplt, argv[++i]);
		}
		else if (strncmp(sz, "-dump_stride", 12) == 0)
		{
			if (sz[12] == '=')
			{
				ops.dumpStride = atoi(sz + 13);
				if (ops.dumpStride < 1)
				{
					fprintf(stderr, "FATAL ERROR: invalid dump stride.\n");
					return false;
				}
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: missing '=' after -dump_stride.\n");
				return false;
			}
		}
		else if (strncmp(sz, "-dump", 5) == 0)
		{
			ops.dumpLevel = FE_DUMP_MAJOR_ITRS;
			if (sz[5] == '=') ops.dumpLevel = atoi(sz + 6);
			if ((ops.dumpLevel < 0) || (ops.dumpLevel > 3))
			{
				fprintf(stderr, "FATAL ERROR: invalid restart level.\n");
				return false;
			}

			if (i<nargs - 1)
			{
				char* szi = argv[i + 1];
				if (szi[0] != '-')
				{
					// assume this is the name of the dump file
					strcpy(ops.szdmp, argv[++i]);
					bdmp = true;
				}
			}
		}
		else if (strcmp(sz, "-o") == 0)
		{
			blog = true;
			strcpy(ops.szlog, argv[++i]);
		}
		else if (strcmp(sz, "-i") == 0)
		{
			++i;
			const char* szext = strrchr(argv[i], '.');
			if (szext == 0)
			{
				// we assume a default extension of .feb if none is provided
				sprintf(ops.szfile, "%s.feb", argv[i]);
			}
			else strcpy(ops.szfile, argv[i]);
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-s") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-s is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "optimize");
			strcpy(ops.szctrl, argv[++i]);
		}
		else if ((strcmp(sz, "-g") == 0) || (strcmp(sz, "-g1") == 0))
		{
			ops.ndebug = 1;
		}
		else if (strcmp(sz, "-g2") == 0)
		{
			ops.ndebug = 2;
		}
		else if (strcmp(sz, "-nosplash") == 0)
		{
			// don't show the welcome message
			ops.bsplash = false;
		}
		else if (strcmp(sz, "-silent") == 0)
		{
			// no output to screen
			ops.bsilent = true;
		}
		else if (strcmp(sz, "-cnf") == 0)	// obsolete: use -config instead
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-config") == 0)
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
		}
		else if (strncmp(sz, "-task", 5) == 0)
		{
			if (sz[5] != '=') { fprintf(stderr, "command line error when parsing task\n"); return false; }
			strcpy(ops.sztask, sz+6);

			if (i<nargs-1)
			{
				char* szi = argv[i+1];
				if (szi[0] != '-')
				{
					// assume this is a control file for the specified task
					strcpy(ops.szctrl, argv[++i]);
				}
			}
		}
		else if (strcmp(sz, "-break") == 0)
		{
			char szbuf[32]={0};
			strcpy(szbuf, argv[++i]);

			add_break_point(szbuf);
		}
		else if (strcmp(sz, "-info")==0)
		{
			FILE* fp = stdout;
			if ((i<nargs-1) && (argv[i+1][0] != '-'))
			{
				fp = fopen(argv[++i], "wt");
				if (fp == 0) fp = stdout;
			}
			fprintf(fp, "compiled on " __DATE__ "\n");

			char* szver = febio::getVersionString();

#ifdef _DEBUG
			fprintf(fp, "FEBio version  = %s (DEBUG)\n", szver);
#else
			fprintf(fp, "FEBio version  = %s\n", szver);
#endif
			if (fp != stdout) fclose(fp);
		}
		else if (strcmp(sz, "-norun") == 0)
		{
			brun = false;
		}

		else if (strcmp(sz, "-import") == 0)
		{
			if ((i < nargs - 1) && (argv[i+1][0] != '-'))
				strcpy(ops.szimp, argv[++i]);
			else
			{
				fprintf(stderr, "FATAL ERROR: insufficient number of arguments for -import.\n");
				return false;
			}
		}
		else if (sz[0] == '-')
		{
			fprintf(stderr, "FATAL ERROR: Invalid command line option.\n");
			return false;
		}
		else
		{
			// if no input file is given yet, we'll assume this is the input file
			if (ops.szfile[0] == 0)
			{
				const char* szext = strrchr(sz, '.');
				if (szext == 0)
				{
					// we assume a default extension of .feb if none is provided
					sprintf(ops.szfile, "%s.feb", sz);
				}
				else
				{
					strcpy(ops.szfile, sz);
				}
				ops.binteractive = false;
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option\n");
				return false;
			}
		}
	}

	// do some sanity checks
	if (strcmp(ops.sztask, "optimize") == 0)
	{
		// make sure we have an input file
		if (ops.szfile[0]==0)
		{
			fprintf(stderr, "FATAL ERROR: no model input file was defined (use -i to define the model input file)\n\n");
			return false;
		}
	}

	// if no task is defined, we assume a std solve is wanted
	if (ops.sztask[0] == 0) strcpy(ops.sztask, "solve");

	// derive the other filenames
	if (ops.szfile[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		char* ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		char szlogbase[256];
		if (ops.szctrl[0])
		{
			strcpy(szlogbase, ops.szctrl);
			ch = strrchr(szlogbase, '.');
			if (ch) *ch = 0;
		}
		else strcpy(szlogbase, szbase);

		if (!blog) sprintf(ops.szlog, "%s.log", szlogbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}
	else if (ops.szctrl[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		strcpy(szbase, ops.szctrl);
		char* ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		if (!blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}

	return brun;
}
