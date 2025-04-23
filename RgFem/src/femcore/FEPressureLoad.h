#pragma once
#include "femcore/FESurfaceLoad.h"
#include  "femcore/FEModelParam.h"

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceLoad
{
public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	//! calculate residual
	void LoadVector(FEGlobalVector& R) override;

	//! calculate stiffness
	void StiffnessMatrix(FELinearSystem& LS) override;

protected:
	FEParamDouble	m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation
	bool			m_blinear;	//!< is the load linear (i.e. it will be calculated in the reference frame and assummed deformation independent)
	bool			m_bshellb;	//!< flag for prescribing pressure on shell bottom

	DECLARE_PARAM_LIST();
};
