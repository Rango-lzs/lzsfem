#include "FESurfacePairConstraint.h"

DEFINE_META_CLASS(FESurfacePairConstraint, FEStepComponent, "");

//-----------------------------------------------------------------------------
FESurfacePairConstraint::FESurfacePairConstraint(FEModel* pfem) : FEStepComponent(pfem)
{
}

//-----------------------------------------------------------------------------
// allocate equations for lagrange multipliers
// (should return the number of equations to be allocated)
int FESurfacePairConstraint::InitEquations(int neq) { return 0; }

//-----------------------------------------------------------------------------
void FESurfacePairConstraint::Update(std::vector<double>& ui) {}
