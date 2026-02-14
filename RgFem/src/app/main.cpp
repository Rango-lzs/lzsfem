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
    // RAII for MPI management to ensure proper cleanup
    struct MPIManager {
        bool initialized = false;
        
        MPIManager(int* argc, char*** argv) {
            int mpi_initialized;
            MPI_Initialized(&mpi_initialized);
            if (!mpi_initialized) {
                int result = MPI_Init(argc, argv);
                if (result == MPI_SUCCESS) {
                    initialized = true;
                }
            }
        }
        
        ~MPIManager() {
            int mpi_finalized;
            MPI_Finalized(&mpi_finalized);
            if (!mpi_finalized && initialized) {
                MPI_Finalize();
            }
        }
    } mpi_manager(&argc, &argv);
#endif

RgFemApp* theApp = RgFemApp::Instance();
    if (!theApp) {
        return 1;
    }

    // initialize the app
    if (!theApp->Init(argc, argv)) {
        return 1;
    }

    // start the app
    int ret = theApp->Run();

    // do some cleanup
    theApp->Finish();

    return ret;
}
