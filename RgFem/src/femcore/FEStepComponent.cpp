#include "FEStepComponent.h"

FEStepComponent::FEStepComponent(FEModel* fem) : FEModelComponent(fem)
{
	// initialize parameters
	m_bactive = true;
}

//-----------------------------------------------------------------------------
bool FEStepComponent::IsActive() const
{
	return m_bactive;
}

//-----------------------------------------------------------------------------
void FEStepComponent::Activate()
{
	m_bactive = true;
}

//-----------------------------------------------------------------------------
void FEStepComponent::Deactivate()
{
	m_bactive = false;
}

//-----------------------------------------------------------------------------
void FEStepComponent::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar& m_bactive;
}
