#include "FESurfacePair.h"
#include "FEFacetSet.h"
#include "FEMesh.h"
#include "basicio/DumpStream.h"

//--------------------------------------------------------
FESurfacePair::FESurfacePair(FEMesh* pm) : m_mesh(pm)
{
	m_surface1 = 0;
	m_surface2 = 0;
}

void FESurfacePair::SetName(const std::string& name)
{
	m_name = name;
}

const std::string& FESurfacePair::GetName() const
{
	return m_name;
}

FEFacetSet* FESurfacePair::GetPrimarySurface()
{
	return m_surface1;
}

void FESurfacePair::SetPrimarySurface(FEFacetSet* pf)
{
	m_surface1 = pf;
}

FEFacetSet* FESurfacePair::GetSecondarySurface()
{
	return m_surface2;
}

void FESurfacePair::SetSecondarySurface(FEFacetSet* pf)
{
	m_surface2 = pf;
}

void FESurfacePair::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	ar& m_surface1;
	ar& m_surface2;
	ar& m_mesh;
}
