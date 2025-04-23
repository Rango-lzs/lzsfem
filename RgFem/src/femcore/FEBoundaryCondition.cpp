#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEFacetSet.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEBoundaryCondition::FEBoundaryCondition(FEModel* pfem) : FEStepComponent(pfem), m_dof(pfem)
{
}

//-----------------------------------------------------------------------------
FEBoundaryCondition::~FEBoundaryCondition()
{
}

//-----------------------------------------------------------------------------
//! fill the prescribed values
void FEBoundaryCondition::PrepStep(std::vector<double>& u, bool brel)
{

}

void FEBoundaryCondition::Serialize(DumpStream& ar)
{
	FEStepComponent::Serialize(ar);
	if (ar.IsShallow() == false) ar & m_dof;
}

//-----------------------------------------------------------------------------
void FEBoundaryCondition::SetDOFList(int ndof)
{
	m_dof.Clear();
	m_dof.AddDof(ndof);
}

//-----------------------------------------------------------------------------
void FEBoundaryCondition::SetDOFList(const std::vector<int>& dofs)
{
	m_dof = dofs;
}

void FEBoundaryCondition::SetDOFList(const FEDofList& dofs)
{
	m_dof = dofs;
}
