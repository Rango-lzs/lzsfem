#include "femcore/FESurfacePairConstraintNL.h"

//-----------------------------------------------------------------------------
FESurfacePairConstraintNL::FESurfacePairConstraintNL(FEModel* pfem) : FENLConstraint(pfem)
{
}

//-----------------------------------------------------------------------------
// allocate equations for lagrange multipliers
// (should return the number of equations to be allocated)
int FESurfacePairConstraintNL::InitEquations(int neq) { return 0; }

//-----------------------------------------------------------------------------
void FESurfacePairConstraintNL::Update(std::vector<double>& ui) {}
