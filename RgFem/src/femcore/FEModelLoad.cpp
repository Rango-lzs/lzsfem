#include "femcore/FEModelLoad.h"
#include "femcore/Solver/FESolver.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEModelLoad::FEModelLoad(FEModel* pfem) : FEStepComponent(pfem), m_dof(pfem)
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
