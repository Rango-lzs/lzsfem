#include "femcore/FENLConstraint.h"

DEFINE_META_CLASS(FENLConstraint, FEStepComponent,"");
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
