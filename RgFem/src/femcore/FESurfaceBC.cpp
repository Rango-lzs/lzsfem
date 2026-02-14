#include "FESurfaceBC.h"
#include "FESurface.h"

DEFINE_META_CLASS(FESurfaceBC, FEBoundaryCondition, "");

FESurfaceBC::FESurfaceBC() : FEBoundaryCondition()
{
	m_surface = nullptr;
}

void FESurfaceBC::SetSurface(FESurface* surface)
{
	m_surface = surface;
}

FESurface* FESurfaceBC::GetSurface()
{
	return m_surface;
}

bool FESurfaceBC::Init()
{
    /*if (m_surface == nullptr) return false;
    if (m_surface->Init() == false) return false;*/
	return FEBoundaryCondition::Init();
}
