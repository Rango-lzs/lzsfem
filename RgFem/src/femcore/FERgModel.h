/*********************************************************************
 * \file   FERgModel.h
 * \brief
 *
 * \author Leizs
 * \date   December 2024
 *********************************************************************/
#pragma once
#include "femcore/fem_export.h"
#include "femcore/FEModel.h"
#include "logger/Logfile.h"
#include "femcore/Timer.h"


//-----------------------------------------------------------------------------
// Dump level determines the times the restart file is written
enum FE_Dump_Level
{
    FE_DUMP_NEVER,       // never write a dump file
    FE_DUMP_MAJOR_ITRS,  // create a dump file at the end of each converged time step
    FE_DUMP_STEP,        // create a dump file at the end of an analysis step
    FE_DUMP_MUST_POINTS  // create a dump file only on must-points
};

//-----------------------------------------------------------------------------
struct ModelStats
{
    int ntimeSteps;     //!< total nr of time steps
    int ntotalIters;    //!< total nr of equilibrium iterations
    int ntotalRHS;      //!< total nr of right hand side evaluations
    int ntotalReforms;  //!< total nr of stiffness reformations
};

class PlotFile;

class FEM_EXPORT FERgModel : public FEModel
{
public:
    FERgModel();
    ~FERgModel();
    bool Init() override;

public:
    //input data from file
    bool Input(const char* szfile);

    void Log(int ntag, const char* szmsg) override;
    // get the log file
    Logfile& GetLogFile()
    {
        return m_log;
    }

public:
    //Return the total timer
    Timer& GetSolveTimer();
    //return number of seconds of time spent in linear solver
    int GetLinearSolverTime();

private:
    Timer m_InputTime;   //!< timer to track time to read model
    Timer m_InitTime;    //!< timer to track model initialization
    Timer m_IOTimer;     //!< timer to track output (include plot, dump, and data)

    PlotFile* m_plot;    //!< the plot file
    bool m_becho;        //!< echo input to logfile
    int m_ndebug;        //!< debug level flag
    bool m_writeMesh;    //!< write a new mesh section

    bool m_bshowErrors;  //!< print warnings and errors
    int m_logLevel;      //!< output level for log file
    int m_dumpLevel;     //!< level or writing restart file
    int m_dumpStride;    //!< write dump file every nth iterations

private:
    // accumulative statistics
    ModelStats m_stats;

protected:                      // file names
    std::string m_sfile_title;  //!< input file title
    std::string m_sfile;        //!< input file name (= path + title)
    std::string m_splot;        //!< plot output file name
    std::string m_slog;         //!< log output file name
    std::string m_sdump;        //!< dump file name
    std::string m_title;        //!< model title

protected:
    bool m_pltAppendOnRestart;
    int m_lastUpdate;

private:
    Logfile m_log;
    DECLARE_PARAM_LIST();
};
