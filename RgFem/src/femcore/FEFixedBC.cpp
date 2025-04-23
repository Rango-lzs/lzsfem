#include "FEFixedBC.h"
#include "FENodeSet.h"
#include "FENode.h"
#include "basicio/DumpStream.h"
#include "logger/log.h"

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FENodalBC(pfem)
{
}

FEFixedBC::FEFixedBC(FEModel* pfem, int dof, FENodeSet* ps) : FENodalBC(pfem)
{
	SetDOFList(dof);
	SetNodeSet(ps);
}

//-----------------------------------------------------------------------------
//! initialization
bool FEFixedBC::Init()
{
	if (GetNodeSet() == nullptr) return false;
	if (GetDofList().Size() == 0)
	{
		feLogError("No degrees of freedom are constrained in %s", GetName().c_str());
		return false;
	}

	return FENodalBC::Init();
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FENodalBC::Activate();

	FENodeSet& nset = *GetNodeSet();
	int n = nset.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i<n; ++i)
	{
		// make sure we only activate open dof's
		FENode& node = *nset.Node(i);
		for (int j=0; j<dofs; ++j)
		{
			int dofj = m_dof[j];
			if (node.get_bc(dofj) == DOF_OPEN)
			{
				node.setDofState(dofj, DOF_FIXED);
				node.set(dofj, 0.0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Deactivate()
{
	FENodalBC::Deactivate();
	FENodeSet* pns = GetNodeSet();
	int N = pns->Size();
	size_t dofs = m_dof.Size();
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *pns->Node(i);

		// set the dof to open
		for (size_t j = 0; j < dofs; ++j)
		{
			node.setDofState(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::CopyFrom(FEBoundaryCondition* bc)
{
	FEFixedBC* fbc = dynamic_cast<FEFixedBC*>(bc);
	m_dof = fbc->GetDofList();
}

//=============================================================================
BEGIN_PARAM_DEFINE(FEFixedDOF, FEFixedBC)
	ADD_PARAMETER(m_dofs, "dofs", 0, "$(dof_list)");
END_PARAM_DEFINE();

FEFixedDOF::FEFixedDOF(FEModel* fem) : FEFixedBC(fem)
{

}

void FEFixedDOF::SetDOFS(const std::vector<int>& dofs)
{
	m_dofs = dofs;
}

bool FEFixedDOF::Init()
{
	SetDOFList(m_dofs);
	return FEFixedBC::Init();
}
