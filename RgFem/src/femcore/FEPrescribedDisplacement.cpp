#include "FEPrescribedDisplacement.h"
#include "units.h"

DEFINE_META_CLASS(FEPrescribedDisplacement, FEPrescribedDOF, "prescribed displacement");

//=======================================================================================
// NOTE: I'm setting FEBoundaryCondition is the base class since I don't want to pull
//       in the parameters of FEPrescribedDOF. 
BEGIN_PARAM_DEFINE(FEPrescribedDisplacement, FEPrescribedDOF)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:displacement)");
	ADD_PARAMETER(m_scale, "value")->setUnits(UNIT_LENGTH)->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
	ADD_PARAMETER(m_brelative, "relative");
END_PARAM_DEFINE();

FEPrescribedDisplacement::FEPrescribedDisplacement() : FEPrescribedDOF()
{
}
