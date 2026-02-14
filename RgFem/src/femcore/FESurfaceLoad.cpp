#include "femcore/FESurfaceLoad.h"
#include "femcore/FEMesh.h"
#include "basicio/DumpStream.h"
#include "femcore/FEModel.h"
#include "femcore/FESurface.h"

DEFINE_META_CLASS(FESurfaceLoad, FEModelLoad, "");

FESurfaceLoad::FESurfaceLoad() : FEModelLoad()
{
	m_psurf = 0;
}

FESurfaceLoad::~FESurfaceLoad(void)
{

}

//! Set the surface to apply the load to
void FESurfaceLoad::SetSurface(FESurface* ps)
{
	m_psurf = ps; 
}

bool FESurfaceLoad::Init()
{
    /*if (m_psurf == 0) return false;
    if (m_psurf->Init() == false) return false;*/
	return FEModelLoad::Init();
}

void FESurfaceLoad::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	if (ar.IsShallow()) return;

	//ar & m_psurf;

	// the mesh manages surfaces for surface loads
	if (m_psurf && ar.IsLoading())
	{
		/*FEMesh* pm = m_psurf->GetMesh();
		pm->AddSurface(m_psurf);*/
	}
}

void FESurfaceLoad::ForceMeshUpdate()
{
	//GetFEModel()->SetMeshUpdateFlag(true);
}
