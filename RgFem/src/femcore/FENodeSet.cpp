#include "FENodeSet.h"
#include "FEMesh.h"
#include "basicio/DumpStream.h"
#include "FEModel.h"

//=============================================================================
// FENodeSet
//-----------------------------------------------------------------------------
FENodeSet::FENodeSet(FEModel* fem) : FEItemList(&fem->GetMesh()), m_Node(&fem->GetMesh())
{
	SetMesh(&fem->GetMesh());
}

//-----------------------------------------------------------------------------
void FENodeSet::Add(int id)
{
	m_Node.Add(id);
}

//-----------------------------------------------------------------------------
void FENodeSet::Add(const std::vector<int>& ns)
{
	m_Node.Add(ns);
}

//-----------------------------------------------------------------------------
void FENodeSet::Add(const FENodeList& ns)
{
	m_Node.Add(ns);
}

//-----------------------------------------------------------------------------
void FENodeSet::Clear()
{
	m_Node.Clear();
}

//-----------------------------------------------------------------------------
FENode* FENodeSet::Node(int i)
{
	FEMesh* mesh = GetMesh();
	return (mesh ? &mesh->Node(m_Node[i]) : nullptr);
}

//-----------------------------------------------------------------------------
const FENode* FENodeSet::Node(int i) const
{
	FEMesh* mesh = GetMesh();
	return (mesh ? &mesh->Node(m_Node[i]) : nullptr);
}

//-----------------------------------------------------------------------------
void FENodeSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Node;
}

void FENodeSet::SaveClass(DumpStream& ar, FENodeSet* p)
{
}

FENodeSet* FENodeSet::LoadClass(DumpStream& ar, FENodeSet* p)
{
	p = new FENodeSet(&ar.GetFEModel());
	return p;
}
