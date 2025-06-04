#include "FELoadCurve.h"

#include "basicio/DumpStream.h"
#include "FEFunction1D.h"
#include "logger/log.h"

BEGIN_PARAM_DEFINE(FELoadCurve, FELoadController)
ADD_PARAMETER(m_int, "interpolate", 0,
              "LINEAR\0STEP\0SMOOTH\0CUBIC SPLINE\0CONTROL POINTS\0APPROXIMATION\0SMOOTH STEP\0");
ADD_PARAMETER(m_ext, "extend", 0, "CONSTANT\0EXTRAPOLATE\0REPEAT\0REPEAT OFFSET\0");
ADD_PARAMETER(m_points, "points");
END_PARAM_DEFINE();

FELoadCurve::FELoadCurve()
    : FELoadController()
{
    m_int = PointCurve::LINEAR;
    m_ext = PointCurve::CONSTANT;
}

FELoadCurve::FELoadCurve(const FELoadCurve& lc)
    : FELoadController(lc)
{
    m_fnc = lc.m_fnc;
}

void FELoadCurve::operator=(const FELoadCurve& lc)
{
    m_fnc = lc.m_fnc;
}

FELoadCurve::~FELoadCurve()
{
}

bool FELoadCurve::Init()
{
    m_fnc.SetInterpolator(m_int);
    m_fnc.SetExtendMode(m_ext);
    m_fnc.SetPoints(m_points);

    // check points
    if (m_fnc.Points() > 1)
    {
        for (int i = 1; i < m_points.size(); ++i)
        {
            double t0 = m_points[i - 1].x();
            double t1 = m_points[i].x();
            if (t0 == t1)
                feLogWarning("Repeated time coordinate in load controller %d", GetID() + 1);
        }
    }

    if (m_fnc.Update() == false)
        return false;
    return FELoadController::Init();
}

void FELoadCurve::Serialize(DumpStream& ar)
{
    FELoadController::Serialize(ar);
    if (ar.IsShallow())
        return;

    if (ar.IsSaving())
    {
        ar << m_int << m_ext;
        ar << m_points;
    }
    else
    {
        ar >> m_int >> m_ext;
        ar >> m_points;

        m_fnc.Clear();
        m_fnc.SetInterpolator(m_int);
        m_fnc.SetExtendMode(m_ext);
        m_fnc.SetPoints(m_points);
        m_fnc.Update();
    }
}

//! evaluates the loadcurve at time
double FELoadCurve::GetValue(double time)
{
    return m_fnc.value(time);
}

bool FELoadCurve::CopyFrom(FELoadCurve* lc)
{
    m_int = lc->m_int;
    m_ext = lc->m_ext;
    m_points = lc->m_points;

    m_fnc = lc->m_fnc;
    return true;
}

void FELoadCurve::Add(double time, double value)
{
    //	m_fnc.Add(time, value);
    m_points.push_back(Vector2d(time, value));
}

void FELoadCurve::Clear()
{
    m_fnc.Clear();
}

void FELoadCurve::SetInterpolation(PointCurve::INTFUNC f)
{
    m_int = f;
    //	m_fnc.SetInterpolator(f);
}

void FELoadCurve::SetExtendMode(PointCurve::EXTMODE f)
{
    m_ext = f;
    //	m_fnc.SetExtendMode(f);
}
