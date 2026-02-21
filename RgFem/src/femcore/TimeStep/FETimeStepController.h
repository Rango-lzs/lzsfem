#pragma once
#include "femcore/FEObjectBase.h"

class AnalysisStep;
class DumpStream;
class FEModel;

//-------------------------------------------------------------------
// Class to control the time step.
//
// Refactored: replaces FEAnalysis* dependency with AnalysisStep* + FEModel*.
// Time parameters (dt, dt0, tend) are now read/written through
// AnalysisStep::getStepControl() / setStepControl() and local state.
//-------------------------------------------------------------------
class FEM_EXPORT FETimeStepController : public FEObjectBase
{
    DECLARE_META_CLASS(FETimeStepController, FEObjectBase);

public:
    FETimeStepController();

    //! Set the analysis step and model
    void SetAnalysis(AnalysisStep* step, FEModel* fem);

    // initialization
    bool Init() override;

    //! reset
    void Reset();

    //! serialize
    void Serialize(DumpStream& ar) override;

    //! copy from
    void CopyFrom(FETimeStepController* tc);

public:
    //! Do a running restart (retry with smaller dt)
    void Retry();

    //! Update Time step based on convergence info
    void AutoTimeStep(int niter);

    //! Adjust for must points
    double CheckMustPoints(double t, double dt);

    //--- Accessors for current time step state ---

    //! Get the current time step size
    double GetTimeStep() const
    {
        return m_dt;
    }

    //! Set the current time step size
    void SetTimeStep(double dt)
    {
        m_dt = dt;
    }

    //! Get the end time for this step
    double GetEndTime() const
    {
        return m_tend;
    }

private:
    AnalysisStep* m_step;  //!< the analysis step (for StepControl access)
    FEModel* m_fem;        //!< the FE model

    // Local copies of time parameters (synced from StepControl)
    double m_dt;    //!< current time step size
    double m_dt0;   //!< initial time step size
    double m_tend;  //!< end time for this step

public:
    int m_nretries;                     //!< nr of retries tried so far
    int m_maxretries;                   //!< max nr of retries allowed per time step
    int m_naggr;                        //!< aggressivness parameter
    int m_nmplc;                        //!< must point load curve number
    int m_nmust;                        //!< current must-point
    int m_next_must;                    //!< next must-point to visit
    int m_iteopt;                       //!< optimum nr of iterations
    double m_dtmin;                     //!< min time step size
    double m_dtmax;                     //!< max time step size
    double m_cutback;                   //!< cut back factor used in aggressive time stepping

    std::vector<double> m_must_points;  //!< the list of must-points
    bool m_mp_repeat;                   //!< repeat must-points
    double m_mp_toff;                   //!< offset for repeat must-points

private:
    double m_ddt;    //!< used by auto-time stepper
    double m_dtp;    //!< previous time step size

    bool m_dtforce;  //!< force max time step

    DECLARE_PARAM_LIST();
};