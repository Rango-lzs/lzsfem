#pragma once
#include "femcore/FEPrescribedBC.h"
#include "femcore/FEModelParam.h"

//-----------------------------------------------------------------------------
//! Boundary condition for prescribing a degree of freedom
class FEM_EXPORT FEPrescribedDOF : public FEPrescribedNodeSet
{
    DECLARE_META_CLASS(FEPrescribedDOF, FEPrescribedNodeSet);

public:
	FEPrescribedDOF();
	FEPrescribedDOF(FEModel* pfem, int dof, FENodeSet* nset);

	void SetDOF(int ndof);
	bool SetDOF(const char* szdof);

	FEPrescribedDOF& SetScale(double s, int lc = -1);

	bool Init() override;

	void CopyFrom(FEBoundaryCondition* pbc) override;

public:
	void GetNodalValues(int n, std::vector<double>& val) override;

protected:
	int				m_dof;		//!< degree of freedom to prescribe
	FEParamDouble	m_scale;	//!< overall scale factor

	DECLARE_PARAM_LIST();
};
