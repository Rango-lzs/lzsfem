#include "femcore/FEModelLoad.h"
#include "femcore/Solver/FESolver.h"
#include "basicio/DumpStream.h"

DEFINE_META_CLASS(FEModelLoad, FEStepComponent, "");

//-----------------------------------------------------------------------------
FEModelLoad::FEModelLoad() : FEStepComponent(), m_dof()
{
}

//-----------------------------------------------------------------------------
const FEDofList& FEModelLoad::GetDofList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEModelLoad::Serialize(DumpStream& ar)
{
	FEStepComponent::Serialize(ar);
	ar & m_dof;
}
