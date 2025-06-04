#include "femcore/FENLConstraint.h"

DEFINE_META_CLASS(FENLConstraint, FEStepComponent, "");
//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint() : FEStepComponent()
{
	static int ncount = 1;
	SetID(ncount++);
}

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint()
{
}
