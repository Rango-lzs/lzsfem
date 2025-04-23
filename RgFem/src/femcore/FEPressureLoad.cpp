#include "femcore/FEPressureLoad.h"
#include "femcore/FEFacetSet.h"
#include "femcore/units.h"
#include "femcore/FEParam.h"
#include "femcore/FESurface.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAM_DEFINE(FEPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_pressure, "pressure")->setUnits(UNIT_PRESSURE)->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
	ADD_PARAMETER(m_bsymm   , "symmetric_stiffness");
	ADD_PARAMETER(m_blinear , "linear");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
END_PARAM_DEFINE()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_pressure = 0.0;
	m_bsymm = true;
	m_bshellb = false;
	m_blinear = false;
}

//-----------------------------------------------------------------------------
bool FEPressureLoad::Init()
{
	// get the degrees of freedom
	m_dof.Clear();
	if (m_bshellb == false)
	{
		//m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
	else
	{
		//m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	}
	if (m_dof.IsEmpty()) return false;

	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FEPressureLoad::LoadVector(FEGlobalVector& R)
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	surf.LoadVector(R, m_dof, m_blinear, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
		
		// evaluate pressure at this material point
		double P = -m_pressure(pt);
		if (m_bshellb) P = -P;

		double J = (pt.dxr ^ pt.dxs).norm();

		// force vector
		Vector3d N = (pt.dxr ^ pt.dxs); N.unit();
		Vector3d t = N*P;

		double H_u = dof_a.shape;

		val[0] = H_u*t.x*J;
		val[1] = H_u*t.y*J;
		val[2] = H_u*t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FEPressureLoad::StiffnessMatrix(FELinearSystem& LS)
{
	// Don't calculate stiffness for a linear load
	if (m_blinear) return;

	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, Matrix& kab) {

		// evaluate pressure at this material point
		double P = -m_pressure(mp);
		if (m_bshellb) P = -P;

		double H_i  = dof_a.shape;
		double Gr_i = dof_a.shape_deriv_r;
		double Gs_i = dof_a.shape_deriv_s;

		double H_j  = dof_b.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		Vector3d vab(0,0,0);
		if (m_bsymm)
			vab = (mp.dxr*(H_j * Gs_i - H_i * Gs_j) - mp.dxs*(H_j * Gr_i - H_i * Gr_j)) * 0.5*P;
		else  
			vab = (mp.dxs*Gr_j - mp.dxr*Gs_j)*(P*H_i);

		Matrix3da K(vab);
		//kab.set(0, 0, K);
	});
}
