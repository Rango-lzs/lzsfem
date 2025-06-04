
#pragma once
#include "femcore/FEObjectBase.h"
class FEAnalysis;
class DumpStream;
class FEModel;

//-------------------------------------------------------------------
// Class to control the time step
class FEM_EXPORT FETimeStepController : public FEObjectBase
{
    DECLARE_META_CLASS(FETimeStepController, FEObjectBase);

public:
    FETimeStepController();

    void SetAnalysis(FEAnalysis* step);

    // initialization
    bool Init() override;

    //! reset
    void Reset();

    //! serialize
    void Serialize(DumpStream& ar) override;

    //! copy from
    void CopyFrom(FETimeStepController* tc);

public:
    //! Do a running restart
    void Retry();

    //! Update Time step
    void AutoTimeStep(int niter);

    //! Adjust for must points
    double CheckMustPoints(double t, double dt);

private:
    FEAnalysis* m_step;

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
