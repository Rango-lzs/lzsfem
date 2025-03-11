#include "FENodeList.h"
#include "FEMesh.h"
#include "DumpStream.h"

FENodeList::FENodeList(FEMesh* mesh) : m_mesh(mesh)
{

}

FENodeList::FENodeList(const FENodeList& nodeList)
{
	m_mesh = nodeList.m_mesh;
	m_nodes = nodeList.m_nodes;
}

FENodeList& FENodeList::operator = (const FENodeList& nodeList)
{
	m_mesh = nodeList.m_mesh;
	m_nodes = nodeList.m_nodes;
	return *this;
}

void FENodeList::Add(int n)
{
	assert(m_mesh);
	assert((n >= 0) && (n < m_mesh->Nodes()));
	m_nodes.push_back(n);
}

void FENodeList::Add(const std::vector<int>& nodeList)
{
	assert(m_mesh);
	m_nodes.insert(m_nodes.end(), nodeList.begin(), nodeList.end());
}

void FENodeList::Add(const FENodeList& nodeList)
{
	assert(m_mesh == nodeList.m_mesh);
	Add(nodeList.m_nodes);
}

void FENodeList::Clear()
{
	m_nodes.clear();
}

FENode* FENodeList::Node(int i)
{
	assert(m_mesh);
	return &m_mesh->Node(m_nodes[i]);
}

const FENode* FENodeList::Node(int i) const
{
	assert(m_mesh);
	return &m_mesh->Node(m_nodes[i]);
}

int FENodeList::Size() const
{
	return (int)m_nodes.size();
}

void FENodeList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false) ar & m_mesh;
	ar & m_nodes;
}

int FENodeList::GlobalToLocalID(int globalId) const
{
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		if (m_nodes[i] == globalId) return i;
	}
	return -1;
}
