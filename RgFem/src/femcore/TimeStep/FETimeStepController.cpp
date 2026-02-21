#include "FETimeStepController.h"

#include "basicio/DumpStream.h"
#include "datastructure/MathUtils.h"
#include "femcore/AnalysisStep/AnalysisStep.h"
#include "femcore/AnalysisStep/StepControl.h"
#include "femcore/FELoadCurve.h"
#include "femcore/FEModel.h"
#include "femcore/FEPointFunction.h"
#include "logger/log.h"

DEFINE_META_CLASS(FETimeStepController, FEObjectBase, "default");

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FETimeStepController, FEParamObject)
ADD_PARAMETER(m_maxretries, "max_retries")->setLongName("max retries");
ADD_PARAMETER(m_iteopt, "opt_iter")->setLongName("optimal iterations");
ADD_PARAMETER(m_dtmin, "dtmin")->setLongName("min stepsize");
ADD_PARAMETER(m_dtmax, "dtmax")->setLongName("max stepsize");
ADD_PARAMETER(m_naggr, "aggressiveness");
ADD_PARAMETER(m_cutback, "cutback");
ADD_PARAMETER(m_dtforce, "dtforce");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FETimeStepController::FETimeStepController()
    : FEObjectBase()
{
    m_step = nullptr;
    m_fem = nullptr;

    m_dt = 0;
    m_dt0 = 0;
    m_tend = 0;

    m_nretries = 0;
    m_maxretries = 5;
    m_naggr = 0;
    m_cutback = 0.5;
    m_nmust = -1;
    m_next_must = -1;
    m_nmplc = -1;
    m_iteopt = 11;
    m_dtmin = 0;
    m_dtmax = 0.1;

    m_ddt = 0;
    m_dtp = 0;
    m_mp_repeat = false;
    m_mp_toff = 0.0;

    m_dtforce = false;
}

//-----------------------------------------------------------------------------
void FETimeStepController::SetAnalysis(AnalysisStep* step, FEModel* fem)
{
    m_step = step;
    m_fem = fem;
}

//-----------------------------------------------------------------------------
//! copy from
void FETimeStepController::CopyFrom(FETimeStepController* tc)
{
    assert(m_step);

    m_naggr = tc->m_naggr;
    m_nmplc = tc->m_nmplc;
    m_iteopt = tc->m_iteopt;
    m_dtmin = tc->m_dtmin;
    m_dtmax = tc->m_dtmax;
    m_cutback = tc->m_cutback;

    m_ddt = tc->m_ddt;
    m_dtp = tc->m_dtp;

    m_must_points = tc->m_must_points;
}

//-----------------------------------------------------------------------------
// initialization
bool FETimeStepController::Init()
{
    // make sure we have a step assigned
    if (m_step == nullptr || m_fem == nullptr)
        return false;

    // Sync local time parameters from AnalysisStep's StepControl
    const StepControl& ctrl = m_step->getStepControl();
    m_dt0 = ctrl.initialTimeIncrement;
    m_dt = ctrl.initialTimeIncrement;
    m_tend = ctrl.totalTime;

    // Also sync dtmin/dtmax from StepControl if not already set by user
    if (m_dtmin <= 0)
        m_dtmin = ctrl.minimumTimeIncrement;
    if (m_dtmax <= 0)
        m_dtmax = ctrl.maximumTimeIncrement;

    // initialize "previous" time step
    m_dtp = m_dt0;

    return true;
}

//-----------------------------------------------------------------------------
//! reset
void FETimeStepController::Reset()
{
    m_dtp = m_dt0;
    m_dt = m_dt0;
    m_nmust = -1;
    m_next_must = -1;
    m_mp_toff = 0.0;
}

//-----------------------------------------------------------------------------
//! Restores data for a running restart

void FETimeStepController::Retry()
{
    feLogEx(m_fem, "Retrying time step. Retry attempt %d of max %d\n\n", m_nretries + 1, m_maxretries);

    // adjust time step
    double dt = m_dt;
    if (m_nretries == 0)
        m_ddt = (dt) / (m_maxretries + 1);

    double dtn;
    if (m_naggr == 0)
        dtn = dt - m_ddt;
    else
        dtn = dt * m_cutback;

    feLogEx(m_fem, "\nAUTO STEPPER: retry step, dt = %lg\n\n", dtn);

    // increase retry counter
    m_nretries++;

    // the new time step cannot be a must-point
    if (m_nmust != -1)
    {
        // if we were at a must-point, make sure we can hit this must-point again
        m_next_must--;
        m_nmust = -1;
    }

    m_dtp = dtn;
    m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! Adjusts the time step size based on the convergence information.
//!	If the previous time step was able to converge in less than
//! m_iteopt iterations the step size is increased, else it is decreased.

void FETimeStepController::AutoTimeStep(int niter)
{
    double dt = m_dt;
    double dtn = m_dtp;
    double told = m_fem->GetCurrentTime();

    // make sure the timestep size is at least the minimum
    if (dtn < m_dtmin)
        dtn = m_dtmin;

    // get the max time step
    double dtmax = m_dtmax;

    // If we have a must-point load curve
    // we take the max step size from the lc
    if (m_nmplc >= 0)
    {
        // TODO: implement must-point load curve support
    }

    // adjust time step size
    if (m_dtforce)
    {
        // if the force flag is set, we just set the time step to the max value
        dtn = dtmax;
    }
    else if (niter > 0)
    {
        double scale = sqrt((double)m_iteopt / (double)niter);

        // Adjust time step size
        if (scale >= 1)
        {
            dtn = dtn + (dtmax - dtn) * MIN(.20, scale - 1);
            dtn = MIN(dtn, 5.0 * m_dtp);
            if (dtmax > 0)
                dtn = MIN(dtn, dtmax);
        }
        else
        {
            dtn = dtn - (dtn - m_dtmin) * (1 - scale);
            if (m_dtmin > 0)
                dtn = MAX(dtn, m_dtmin);
            if (dtmax > 0)
                dtn = MIN(dtn, dtmax);
        }
    }
    else if (niter == 0)
    {
        if (m_dtmin > 0)
            dtn = MAX(dtn, m_dtmin);
        if (dtmax > 0)
            dtn = MIN(dtn, dtmax);
    }

    // Report new time step size
    if (dtn > dt)
        feLogEx(m_fem, "\nAUTO STEPPER: increasing time step, dt = %lg\n\n", dtn);
    else if (dtn < dt)
        feLogEx(m_fem, "\nAUTO STEPPER: decreasing time step, dt = %lg\n\n", dtn);

    // Store this time step value. This is the value that will be used to evaluate
    // the next time step increment. This will not include adjustments due to the must-point
    // controller since this could create really small time steps that may be difficult to
    // recover from.
    m_dtp = dtn;

    // check for mustpoints
    if (m_must_points.empty() == false)
        dtn = CheckMustPoints(told, dtn);

    // make sure we are not exceeding the final time
    if (told + dtn > m_tend)
    {
        dtn = m_tend - told;
        feLogEx(m_fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtn);
    }

    // store time step size
    assert(dtn > 0);
    m_dt = dtn;
}

//-----------------------------------------------------------------------------
//! This function makes sure that no must points are passed. It returns an
//! updated value (less than dt) if t + dt would pass a must point. Otherwise
//! it returns dt.
//! \param t current time
//! \param dt current time step
//! \return updated time step.
double FETimeStepController::CheckMustPoints(double t, double dt)
{
    const double eps = m_tend * 1e-12;

    m_nmust = -1;
    const int points = (int)m_must_points.size();

    if (m_next_must >= points)
    {
        if (m_mp_repeat)
        {
            m_mp_toff += m_must_points.back();
            m_next_must = -1;
        }
        else
            return dt;
    }

    // set the first must-point if it has not been set
    if (m_next_must < 0)
    {
        m_next_must = 0;
        while ((m_next_must < points) && (m_must_points[m_next_must] + m_mp_toff < t + eps))
        {
            m_next_must++;
            if (m_next_must >= points)
            {
                if (m_mp_repeat)
                {
                    m_next_must = 0;
                    m_mp_toff += m_must_points.back();
                }
                else
                    return dt;
            }
        }
    }

    assert(m_next_must < points);
    double tmust = m_must_points[m_next_must] + m_mp_toff;
    assert(tmust + eps > t);

    double dtnew = dt;
    double tnew = t + dt;
    if (tmust < tnew + eps)
    {
        dtnew = tmust - t;
        feLogEx(m_fem, "MUST POINT CONTROLLER: adjusting time step. dt = %lg\n\n", dtnew);
        m_nmust = m_next_must++;
    }

    return dtnew;
}

//-----------------------------------------------------------------------------
//! serialize
void FETimeStepController::Serialize(DumpStream& ar)
{
    FEObjectBase::Serialize(ar);
    ar& m_nretries;
    ar& m_nmplc;
    ar& m_nmust;
    ar& m_next_must;
    ar& m_mp_toff;
    ar& m_mp_repeat;
    ar& m_ddt& m_dtp;
    ar& m_dt& m_dt0& m_tend;
    ar& m_must_points;
}