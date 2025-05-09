#pragma once

#include "FEBoundaryCondition.h"

class FESurface;

class FEM_EXPORT FESurfaceBC : public FEBoundaryCondition
{
    DECLARE_META_CLASS(FESurfaceBC, FEBoundaryCondition);

public:
	FESurfaceBC(FEModel* fem);

	bool Init() override;

	void SetSurface(FESurface* surface);

	FESurface* GetSurface();

private:
	FESurface* m_surface;
};
