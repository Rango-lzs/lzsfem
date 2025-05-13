#include "FERgModel.h"

#include "input/febio/FEBioImport.h"
#include "input/febio/FEBioModelBuilder.h"
#include "femcore/Callback.h"
#include "logger/log.h"

#include <fstream>
#include <iostream>
#include <sstream>

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FERgModel, FEModel)
ADD_PARAMETER(m_title, "title");
ADD_PARAMETER(m_logLevel, "log_level");
END_PARAM_DEFINE();


// Constructor of FERgModel class.
FERgModel::FERgModel()
{
    m_logLevel = 1;

    m_dumpLevel = FE_DUMP_NEVER;
    m_dumpStride = 1;

    // --- I/O-Data ---
    m_ndebug = 0;
    m_becho = true;
    m_plot = nullptr;
    m_writeMesh = false;

    m_stats.ntimeSteps = 0;
    m_stats.ntotalIters = 0;
    m_stats.ntotalRHS = 0;
    m_stats.ntotalReforms = 0;

    m_pltAppendOnRestart = true;

    m_lastUpdate = -1;

    m_bshowErrors = true;

    // Add the output callback
    // We call this function always since we want to flush the logfile for each event.
    // 
    //Rango TODO:
    //AddCallback(handleCB, CB_ALWAYS, this);
}

//-----------------------------------------------------------------------------
FERgModel::~FERgModel()
{
    // close the plot file
    if (m_plot)
    {
        delete m_plot;
        m_plot = 0;
    }
    m_log.close();
}



//-----------------------------------------------------------------------------
Timer& FERgModel::GetSolveTimer()
{
    return *GetTimer(Timer_ModelSolve);
}

//-----------------------------------------------------------------------------
//! return number of seconds of time spent in linear solver
int FERgModel::GetLinearSolverTime()
{
    Timer* t = GetTimer(TimerID::Timer_LinSolve);
    return (int)t->peek();
}

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FERgModel::Init()
{
    TimerTracker t(&m_InitTime);

    //// Open the logfile
    //if (m_logLevel != 0)
    //{
    //    if (InitLogFile() == false)
    //        return false;
    //}

    //FEBioPlotFile* pplt = nullptr;
    //m_lastUpdate = -1;

    //// open plot database file
    //FEAnalysis* step = GetCurrentStep();
    //if (step->GetPlotLevel() != FE_PLOT_NEVER)
    //{
    //    if (m_plot == 0)
    //        InitPlotFile();
    //}

    //// see if a valid dump file name is defined.
    //const std::string& sdmp = GetDumpFileName();
    //if (sdmp.empty())
    //{
    //    // if not, we take the input file name and set the extension to .dmp
    //    char sz[1024] = {0};
    //    strcpy(sz, GetInputFileName().c_str());
    //    char* ch = strrchr(sz, '.');
    //    if (ch)
    //        *ch = 0;
    //    strcat(sz, ".dmp");
    //    SetDumpFilename(sz);
    //}

    //// initialize data records
    //DataStore& dataStore = GetDataStore();
    //for (int i = 0; i < dataStore.Size(); ++i)
    //{
    //    if (dataStore.GetDataRecord(i)->Initialize() == false)
    //        return false;
    //}

    // initialize model data
    if (FEModel::Init() == false)
    {
        feLogError("Model initialization failed");
        return false;
    }

    // Alright, all initialization is done, so let's get busy !
    return true;
}

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FERgModel::Input(const char* szfile)
{
    // start the timer
    TimerTracker t(&m_InputTime);

    // create file reader
    FEBioImport fim;

    // override the default model builder
    fim.SetModelBuilder(new FEBioModelBuilder(*this));

    feLog("Reading file %s ...", szfile);

    // Load the file
    if (fim.Load(*this, szfile) == false)
    {
        feLog("FAILED!\n");
        char szerr[256];
        fim.GetErrorMessage(szerr);
        feLogError(szerr);

        return false;
    }
    else
        feLog("SUCCESS!\n");

    //// set the input file name
    //SetInputFilename(szfile);

    //// see if user redefined output filenames
    //if (fim.m_szdmp[0])
    //    SetDumpFilename(fim.m_szdmp);
    //if (fim.m_szlog[0])
    //    SetLogFilename(fim.m_szlog);
    //if (fim.m_szplt[0])
    //    SetPlotFilename(fim.m_szplt);

    //// add the data records
    //int ND = (int)fim.m_data.size();
    //for (int i = 0; i < ND; ++i)
    //    AddDataRecord(fim.m_data[i]);

    // we're done reading
    return true;
}