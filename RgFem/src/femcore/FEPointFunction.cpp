#include "FEPointFunction.h"
#include "basicio/DumpStream.h"
#include "logger/log.h"
#include "BSpline.h"

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEPointFunction, FEFunction1D)
	ADD_PARAMETER(m_int, "interpolate", 0, "linear\0step\0smooth\0cubic spline\0control points\0approximation\0");
	ADD_PARAMETER(m_ext, "extend"     , 0, "constant\0extrapolate\0repeat\0repeat offset\0");
    ADD_PARAMETER(m_bln, "log")->SetFlags(FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_points, "points");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
//! default constructor
FEPointFunction::FEPointFunction(FEModel* fem) : FEFunction1D(fem)
{
	m_int = PointCurve::LINEAR;
	m_ext = PointCurve::CONSTANT;
    m_bln = false;
}

//-----------------------------------------------------------------------------
FEPointFunction::~FEPointFunction()
{

}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
bool FEPointFunction::Init()
{
	m_fnc.SetInterpolator(m_int);
	m_fnc.SetExtendMode(m_ext);
	m_fnc.SetPoints(m_points);
	if (m_fnc.Update() == false) return false;

    return FEFunction1D::Init();
}

//-----------------------------------------------------------------------------
//! Clears the loadcurve data
void FEPointFunction::Clear()
{ 
	m_fnc.Clear();
}

//-----------------------------------------------------------------------------
//! return nr of points
int FEPointFunction::Points() const
{ 
	return (int) m_points.size(); 
}

//-----------------------------------------------------------------------------
//! set the points
void FEPointFunction::SetPoints(const std::vector<Vector2d>& pts)
{
	m_points = pts;
}

//-----------------------------------------------------------------------------
// Sets the time and data value of point i
// This function assumes that the load curve data has already been created
//
void FEPointFunction::SetPoint(int i, double x, double y)
{
	m_fnc.SetPoint(i, x, y);
}

//-----------------------------------------------------------------------------
//! Set the type of interpolation
void FEPointFunction::SetInterpolation(int fnc) { m_int = fnc; }

//-----------------------------------------------------------------------------
//! Set the extend mode
void FEPointFunction::SetExtendMode(int mode) { m_ext = mode; }

//-----------------------------------------------------------------------------
//! returns point i
LOADPOINT FEPointFunction::LoadPoint(int i) const
{ 
	const Vector2d& p = m_points[i];
	LOADPOINT lp;
	lp.time  = p.x();
	lp.value = p.y();
	return lp; 
}

//-----------------------------------------------------------------------------
//! This function adds a datapoint to the loadcurve. The datapoint is inserted
//! at the appropriate place by examining the time parameter.

void FEPointFunction::Add(double x, double y)
{
	// find the place to insert the data point
	int n = 0;
	int nsize = m_points.size();
	while ((n<nsize) && (m_points[n].x() < x)) ++n;

	// insert loadpoint
	m_points.insert(m_points.begin() + n, Vector2d(x, y));
}

//-----------------------------------------------------------------------------
void FEPointFunction::Scale(double s)
{
	for (int i = 0; i < Points(); ++i)
	{
		m_points[i].y() *= s;
	}
}

//-----------------------------------------------------------------------------
double FEPointFunction::value(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.value(time);
}

//-----------------------------------------------------------------------------
void FEPointFunction::Serialize(DumpStream& ar)
{
	FEFunction1D::Serialize(ar);
	if (ar.IsShallow()) return;

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

//-----------------------------------------------------------------------------
double FEPointFunction::derive(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.derive(time);
}

//-----------------------------------------------------------------------------
double FEPointFunction::deriv2(double time) const
{
    if (m_bln) time = (time > 0) ? log(time) : m_points[0].x();
	return m_fnc.deriv2(time);
}

double FEPointFunction::integrate(double a, double b) const
{
	return m_fnc.integrate(a, b);
}

//-----------------------------------------------------------------------------
FEFunction1D* FEPointFunction::copy()
{
	FEPointFunction* f = new FEPointFunction(GetFEModel());

	f->m_int = m_int;
	f->m_ext = m_ext;
	f->m_points = m_points;
	f->m_fnc = m_fnc;
	return f;
}

//-----------------------------------------------------------------------------
void FEPointFunction::CopyFrom(const FEPointFunction& f)
{
	m_int = f.m_int;
	m_ext = f.m_ext;
	m_points = f.m_points;
	m_fnc = f.m_fnc;
}

//-----------------------------------------------------------------------------
void FEPointFunction::CopyFrom(const PointCurve& f)
{
	m_int = f.GetInterpolator();
	m_ext = f.GetExtendMode();
	m_points = f.GetPoints();
	m_fnc = f;
}
