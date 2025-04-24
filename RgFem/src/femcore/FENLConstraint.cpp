#include "femcore/FENLConstraint.h"

//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint(FEModel* pfem) : FEStepComponent(pfem)
{
	static int ncount = 1;
	SetID(ncount++);
}

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint()
{
}
