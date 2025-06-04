#include "FEBodyForce.h"
#include "materials/FESolidMaterial.h"
#include "femcore/Domain/FEElasticDomain.h"

//-----------------------------------------------------------------------------
FEBodyForce::FEBodyForce() : FEBodyLoad()
{
}

//-----------------------------------------------------------------------------
// NOTE: Work in progress! Working on integrating body loads as model loads
void FEBodyForce::LoadVector(FEGlobalVector& R)
{
	for (int i = 0; i<Domains(); ++i)
	{
		FEDomain* dom = Domain(i);
		FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(dom);
		if (edom) edom->BodyForce(R, *this);
	}
}

//-----------------------------------------------------------------------------
// NOTE: Work in progress! Working on integrating body loads as model loads
void FEBodyForce::StiffnessMatrix(FELinearSystem& LS)
{
	for (int i = 0; i<Domains(); ++i)
	{
		FEDomain* dom = Domain(i);
		FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(dom);
		if (edom) edom->BodyForceStiffness(LS, *this);
	}
}
