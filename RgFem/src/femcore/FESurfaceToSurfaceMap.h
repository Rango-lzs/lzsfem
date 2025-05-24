#pragma once
//#include "FEDataGenerator.h"
#include "FEFunction1D.h"
#include "FEClosestPointProjection.h"
#include "FENormalProjection.h"

class FEModel;
class FESurface;

class FESurfaceToSurfaceMap : public FEElemDataGenerator
{
public:
	FESurfaceToSurfaceMap(FEModel* fem);
	~FESurfaceToSurfaceMap();

	bool Init() override;

	void value(const Vector3d& x, double& data) override;

	FEDomainMap* Generate() override;

private:
	FEFunction1D*	m_func;
	
private:
	FESurface*	m_surf1;
	FESurface*	m_surf2;
	FEClosestPointProjection*	m_ccp1;
	FEClosestPointProjection*	m_ccp2;
	bool			m_binverted;

	DECLARE_FECORE_CLASS();
};
