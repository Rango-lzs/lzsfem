#pragma once
#include "FEModelLoad.h"
#include  "datastructure/Vector3d.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FESurface;
class FESurfaceMaterialPoint;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FEM_EXPORT FESurfaceLoad : public FEModelLoad
{
public:
	FESurfaceLoad(FEModel* pfem);
	virtual ~FESurfaceLoad(void);

	//! Set the surface to apply the load to
	virtual void SetSurface(FESurface* ps);

	bool Init() override;

	//! Get the surface
	FESurface& GetSurface() { return *m_psurf; }

	void Serialize(DumpStream& ar) override;

public:
    virtual double ScalarLoad(FESurfaceMaterialPoint& mp) { return 0; }
    virtual Vector3d VectorLoad(FESurfaceMaterialPoint& mp) { return Vector3d(0,0,0); }

	// TODO: Can I get rid of this?
	// This is needed to update the mesh after some surface loads, which aren't really
	// surface loads modify boundary conditions
	void ForceMeshUpdate();

protected:
	FESurface*	m_psurf;

	FECORE_BASE_CLASS(FESurfaceLoad)
};
