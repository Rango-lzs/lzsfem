#include "FEInitialCondition.h"
#include "basicio/DumpStream.h"
#include "materials/FEMaterialPoint.h"
#include "FENode.h"
#include "FEModel.h"

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEStepComponent(pfem)
{
}

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FENodalIC, FEInitialCondition)
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FENodalIC::FENodalIC(FEModel* fem) : FEInitialCondition(fem), m_dofs(fem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
// set the nodeset for this component
void FENodalIC::SetNodeSet(FENodeSet* nset)
{
	m_nodeSet = nset;
}

//-----------------------------------------------------------------------------
// get the node set
FENodeSet* FENodalIC::GetNodeSet()
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
// set the list of degrees of freedom
void FENodalIC::SetDOFList(const FEDofList& dofList)
{
	m_dofs = dofList;
}

//-----------------------------------------------------------------------------
bool FENodalIC::Init()
{
	if (m_nodeSet == nullptr) return false;
	return FEInitialCondition::Init();
}

//-----------------------------------------------------------------------------
void FENodalIC::Activate()
{
	FEStepComponent::Activate();
	if (m_dofs.IsEmpty()) return;

	int dofs = (int)m_dofs.Size();
	std::vector<double> val(dofs, 0.0);

	int N = (int)m_nodeSet->Size();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = *m_nodeSet->Node(i);

		// get the nodal values
		GetNodalValues(i, val);
		
		for (int j = 0; j < dofs; ++j)
		{
			node.set(m_dofs[j], val[j]);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FENodalIC::Serialize(DumpStream& ar)
{
	FEStepComponent::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_dofs;
	ar & m_nodeSet;
}

//======================================================================================
BEGIN_PARAM_DEFINE(FEInitialDOF, FENodalIC)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(m_data, "value");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* pfem) : FENodalIC(pfem)
{
	m_dof = -1;
	m_data = 0.0;
}

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* fem, int ndof, FENodeSet* nset) : FENodalIC(fem)
{
	SetDOF(ndof);
	SetNodeSet(nset);
	m_data = 0.0;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetDOF(int ndof) { m_dof = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialDOF::SetDOF(const char* szdof)
{
	FEModel* fem = GetFEModel();
	int ndof = fem->GetDOFIndex(szdof);
	assert(ndof >= 0);
	if (ndof < 0) return false;
	SetDOF(ndof);
	return true;
}

//-----------------------------------------------------------------------------
bool FEInitialDOF::Init()
{
	if (FENodalIC::Init() == false) return false;
	if (m_dof == -1) return false;
	FEDofList dofs(GetFEModel());
	if (dofs.AddDof(m_dof) == false) return false;
	SetDOFList(dofs);
	return true;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetValue(double v)
{
	m_data = v;
}

//-----------------------------------------------------------------------------
// return the values for node i
void FEInitialDOF::GetNodalValues(int inode, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[inode];
	const FENode& node = *nset.Node(inode);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = inode;

	val[0] = m_data(mp);
}
