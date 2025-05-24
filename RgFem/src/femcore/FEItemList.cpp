#include "FEItemList.h"
#include "basicio/DumpStream.h"
#include "FEModel.h"

FEItemList::FEItemList(FEMesh* mesh)// : FEObjectBase(fem)
{
	m_mesh = mesh;
}

FEItemList::~FEItemList() {}

void FEItemList::Serialize(DumpStream& ar)
{
	ar & m_name;
}

// get the mesh
FEMesh* FEItemList::GetMesh() const
{
	return m_mesh;
}

void FEItemList::SetMesh(FEMesh* mesh)
{
	m_mesh = mesh;
}

const std::string& FEItemList::GetName() const
{
	return m_name;
}

void FEItemList::SetName(const std::string& name)
{
	m_name = name;
}

FEItemList* FEItemList::LoadClass(DumpStream& ar, FEItemList* p)
{
	assert(false);
	return nullptr;
}

void FEItemList::SaveClass(DumpStream& ar, FEItemList* p)
{
	assert(false);
}
