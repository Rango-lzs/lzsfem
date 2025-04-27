/*****************************************************************//**
 * \file   main.cpp
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#include "app/RgFemApp.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

//-----------------------------------------------------------------------------
// The starting point of the application
//
int main(int argc, char* argv[])
{
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif

	RgFemApp* theApp = RgFemApp::Instance();
	if (!theApp)
	{
		return 1;
	}

	// initialize the app
	if (!theApp->Init(argc, argv))
	{
		return 1;
	}

	// start the app
	int ret = theApp->Run();

	// do some cleanup
	theApp->Finalize();

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return ret;
}
