#include "RgLoadController.h"
#include "logger/log.h"
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//=============================================================================
// RgLoadController
//=============================================================================

DEFINE_META_CLASS(RgLoadController, FEObjectBase, "");

RgLoadController::RgLoadController()
    : m_value(0.0), m_type(LC_LOADCURVE), 
      m_interpolation(INTERP_LINEAR), m_extendMode(EXTEND_CONSTANT)
{
}

RgLoadController::~RgLoadController()
{
}

bool RgLoadController::Init()
{
    return FEObjectBase::Init();
}

void RgLoadController::Evaluate(double time)
{
    // Base class does nothing, derived classes override
    m_value = 0.0;
}

void RgLoadController::Serialize(DumpStream& ar)
{
    FEObjectBase::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_value;
    ar & m_type;
    ar & m_interpolation;
    ar & m_extendMode;
}

//=============================================================================
// RgLoadCurve
//=============================================================================

DEFINE_META_CLASS(RgLoadCurve, RgLoadController, "");

RgLoadCurve::RgLoadCurve()
{
    m_type = LC_LOADCURVE;
}

bool RgLoadCurve::Init()
{
    if (!RgLoadController::Init()) return false;
    
    if (m_points.empty())
    {
        feLogWarning("Load curve has no points");
        return false;
    }
    
    // Sort points by time
    std::sort(m_points.begin(), m_points.end(), 
              [](const Point& a, const Point& b) { return a.time < b.time; });
    
    return true;
}

void RgLoadCurve::Evaluate(double time)
{
    if (m_points.empty())
    {
        m_value = 0.0;
        return;
    }
    
    m_value = Interpolate(time);
}

void RgLoadCurve::AddPoint(double time, double value)
{
    m_points.push_back(Point(time, value));
}

void RgLoadCurve::Clear()
{
    m_points.clear();
}

void RgLoadCurve::GetPoint(int i, double& time, double& value) const
{
    if (i >= 0 && i < (int)m_points.size())
    {
        time = m_points[i].time;
        value = m_points[i].value;
    }
}

void RgLoadCurve::SetPoint(int i, double time, double value)
{
    if (i >= 0 && i < (int)m_points.size())
    {
        m_points[i].time = time;
        m_points[i].value = value;
    }
}

void RgLoadCurve::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    if (ar.IsSaving())
    {
        int npoints = (int)m_points.size();
        ar << npoints;
        for (int i = 0; i < npoints; ++i)
        {
            ar << m_points[i].time;
            ar << m_points[i].value;
        }
    }
    else
    {
        int npoints;
        ar >> npoints;
        m_points.resize(npoints);
        for (int i = 0; i < npoints; ++i)
        {
            ar >> m_points[i].time;
            ar >> m_points[i].value;
        }
    }
}

double RgLoadCurve::Interpolate(double t)
{
    if (m_points.empty()) return 0.0;
    if (m_points.size() == 1) return m_points[0].value;
    
    // Check if outside range
    if (t < m_points.front().time || t > m_points.back().time)
    {
        return ExtendValue(t);
    }
    
    // Interpolate based on type
    switch (m_interpolation)
    {
        case INTERP_LINEAR:
            return LinearInterpolate(t);
        case INTERP_STEP:
            return StepInterpolate(t);
        case INTERP_SMOOTH:
            return SmoothInterpolate(t);
        case INTERP_AKIMA:
            // Fall back to smooth for now
            return SmoothInterpolate(t);
        default:
            return LinearInterpolate(t);
    }
}

double RgLoadCurve::LinearInterpolate(double t)
{
    // Find the two points that bracket t
    int n = (int)m_points.size();
    
    for (int i = 0; i < n - 1; ++i)
    {
        if (t >= m_points[i].time && t <= m_points[i + 1].time)
        {
            double t0 = m_points[i].time;
            double t1 = m_points[i + 1].time;
            double v0 = m_points[i].value;
            double v1 = m_points[i + 1].value;
            
            if (fabs(t1 - t0) < 1e-12)
                return v0;
            
            double s = (t - t0) / (t1 - t0);
            return v0 + s * (v1 - v0);
        }
    }
    
    return m_points.back().value;
}

double RgLoadCurve::StepInterpolate(double t)
{
    // Piecewise constant - return value of previous point
    int n = (int)m_points.size();
    
    for (int i = n - 1; i >= 0; --i)
    {
        if (t >= m_points[i].time)
            return m_points[i].value;
    }
    
    return m_points[0].value;
}

double RgLoadCurve::SmoothInterpolate(double t)
{
    // Cubic Hermite interpolation
    int n = (int)m_points.size();
    
    for (int i = 0; i < n - 1; ++i)
    {
        if (t >= m_points[i].time && t <= m_points[i + 1].time)
        {
            double t0 = m_points[i].time;
            double t1 = m_points[i + 1].time;
            double v0 = m_points[i].value;
            double v1 = m_points[i + 1].value;
            
            if (fabs(t1 - t0) < 1e-12)
                return v0;
            
            double s = (t - t0) / (t1 - t0);
            
            // Compute tangents (finite difference)
            double m0, m1;
            
            if (i == 0)
                m0 = (v1 - v0) / (t1 - t0);
            else
                m0 = ((v1 - v0) / (t1 - t0) + 
                      (v0 - m_points[i - 1].value) / (t0 - m_points[i - 1].time)) * 0.5;
            
            if (i == n - 2)
                m1 = (v1 - v0) / (t1 - t0);
            else
                m1 = ((v1 - v0) / (t1 - t0) + 
                      (m_points[i + 2].value - v1) / (m_points[i + 2].time - t1)) * 0.5;
            
            // Hermite basis functions
            double h00 = (1 + 2 * s) * (1 - s) * (1 - s);
            double h10 = s * (1 - s) * (1 - s);
            double h01 = s * s * (3 - 2 * s);
            double h11 = s * s * (s - 1);
            
            return h00 * v0 + h10 * (t1 - t0) * m0 + h01 * v1 + h11 * (t1 - t0) * m1;
        }
    }
    
    return m_points.back().value;
}

double RgLoadCurve::ExtendValue(double t)
{
    if (m_points.empty()) return 0.0;
    
    switch (m_extendMode)
    {
        case EXTEND_CONSTANT:
            // Return end values
            if (t < m_points.front().time)
                return m_points.front().value;
            else
                return m_points.back().value;
            
        case EXTEND_EXTRAPOLATE:
        {
            // Linear extrapolation
            if (t < m_points.front().time)
            {
                if (m_points.size() < 2)
                    return m_points.front().value;
                
                double t0 = m_points[0].time;
                double t1 = m_points[1].time;
                double v0 = m_points[0].value;
                double v1 = m_points[1].value;
                
                double slope = (v1 - v0) / (t1 - t0);
                return v0 + slope * (t - t0);
            }
            else
            {
                int n = (int)m_points.size();
                if (n < 2)
                    return m_points.back().value;
                
                double t0 = m_points[n - 2].time;
                double t1 = m_points[n - 1].time;
                double v0 = m_points[n - 2].value;
                double v1 = m_points[n - 1].value;
                
                double slope = (v1 - v0) / (t1 - t0);
                return v1 + slope * (t - t1);
            }
        }
        
        case EXTEND_REPEAT:
        {
            // Shift time to be within range
            double t0 = m_points.front().time;
            double t1 = m_points.back().time;
            double period = t1 - t0;
            
            if (period < 1e-12)
                return m_points.front().value;
            
            double t_shifted = fmod(t - t0, period);
            if (t_shifted < 0) t_shifted += period;
            t_shifted += t0;
            
            return Interpolate(t_shifted);
        }
        
        case EXTEND_CYCLE:
        {
            // Similar to repeat but cycles the values
            double t0 = m_points.front().time;
            double t1 = m_points.back().time;
            double period = t1 - t0;
            
            if (period < 1e-12)
                return m_points.front().value;
            
            double t_shifted = fmod(t - t0, period);
            if (t_shifted < 0) t_shifted += period;
            t_shifted += t0;
            
            return Interpolate(t_shifted);
        }
        
        default:
            return m_points.front().value;
    }
}

//=============================================================================
// RgLinearController
//=============================================================================

DEFINE_META_CLASS(RgLinearController, RgLoadController, "");

RgLinearController::RgLinearController()
    : m_t0(0.0), m_t1(1.0), m_v0(0.0), m_v1(1.0)
{
    m_type = LC_LINEAR;
}

void RgLinearController::SetParameters(double t0, double t1, double v0, double v1)
{
    m_t0 = t0;
    m_t1 = t1;
    m_v0 = v0;
    m_v1 = v1;
}

void RgLinearController::Evaluate(double time)
{
    if (time <= m_t0)
        m_value = m_v0;
    else if (time >= m_t1)
        m_value = m_v1;
    else
    {
        double s = (time - m_t0) / (m_t1 - m_t0);
        m_value = m_v0 + s * (m_v1 - m_v0);
    }
}

void RgLinearController::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    ar & m_t0 & m_t1 & m_v0 & m_v1;
}

//=============================================================================
// RgStepController
//=============================================================================

DEFINE_META_CLASS(RgStepController, RgLoadController, "");

RgStepController::RgStepController()
    : m_t_step(0.0), m_v0(0.0), m_v1(1.0)
{
    m_type = LC_STEP;
}

void RgStepController::SetParameters(double t_step, double v0, double v1)
{
    m_t_step = t_step;
    m_v0 = v0;
    m_v1 = v1;
}

void RgStepController::Evaluate(double time)
{
    m_value = (time < m_t_step) ? m_v0 : m_v1;
}

void RgStepController::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    ar & m_t_step & m_v0 & m_v1;
}

//=============================================================================
// RgSmoothStepController
//=============================================================================

DEFINE_META_CLASS(RgSmoothStepController, RgLoadController, "");

RgSmoothStepController::RgSmoothStepController()
    : m_t0(0.0), m_t1(1.0), m_v0(0.0), m_v1(1.0)
{
    m_type = LC_SMOOTH_STEP;
}

void RgSmoothStepController::SetParameters(double t0, double t1, double v0, double v1)
{
    m_t0 = t0;
    m_t1 = t1;
    m_v0 = v0;
    m_v1 = v1;
}

void RgSmoothStepController::Evaluate(double time)
{
    if (time <= m_t0)
        m_value = m_v0;
    else if (time >= m_t1)
        m_value = m_v1;
    else
    {
        double s = (time - m_t0) / (m_t1 - m_t0);
        // Smoothstep function: 3s^2 - 2s^3
        double smooth_s = s * s * (3.0 - 2.0 * s);
        m_value = m_v0 + smooth_s * (m_v1 - m_v0);
    }
}

void RgSmoothStepController::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    ar & m_t0 & m_t1 & m_v0 & m_v1;
}

//=============================================================================
// RgSineController
//=============================================================================

DEFINE_META_CLASS(RgSineController, RgLoadController, "");

RgSineController::RgSineController()
    : m_amplitude(1.0), m_frequency(1.0), m_phase(0.0), m_offset(0.0)
{
    m_type = LC_SINE;
}

void RgSineController::SetParameters(double amplitude, double frequency, 
                                      double phase, double offset)
{
    m_amplitude = amplitude;
    m_frequency = frequency;
    m_phase = phase;
    m_offset = offset;
}

void RgSineController::Evaluate(double time)
{
    double omega = 2.0 * M_PI * m_frequency;
    m_value = m_offset + m_amplitude * sin(omega * time + m_phase);
}

void RgSineController::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    ar & m_amplitude & m_frequency & m_phase & m_offset;
}

//=============================================================================
// RgConstantController
//=============================================================================

DEFINE_META_CLASS(RgConstantController, RgLoadController, "");

RgConstantController::RgConstantController()
    : m_constantValue(1.0)
{
    m_type = LC_CONSTANT;
}

void RgConstantController::SetValue(double value)
{
    m_constantValue = value;
    m_value = value;
}

void RgConstantController::Evaluate(double time)
{
    m_value = m_constantValue;
}

void RgConstantController::Serialize(DumpStream& ar)
{
    RgLoadController::Serialize(ar);
    ar & m_constantValue;
}
