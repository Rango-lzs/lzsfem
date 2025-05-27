#include "femcore/FENodalBC.h"
#include "basicio/DumpStream.h"

DEFINE_META_CLASS(FENodalBC, FEBoundaryCondition, "");

//==================================================================
FENodalBC::FENodalBC(FEModel* fem) : FEBoundaryCondition(fem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
void FENodalBC::SetNodeSet(FENodeSet* nodeSet)
{
	m_nodeSet = nodeSet;
}

//-----------------------------------------------------------------------------
FENodeSet* FENodalBC::GetNodeSet()
{
	return m_nodeSet;
}

void FENodalBC::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_nodeSet;
}
