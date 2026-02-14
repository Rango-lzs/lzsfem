#include "FEBodyForce.h"
#include "materials/FESolidMaterial.h"
#include "femcore/Domain/RgSolidDomain.h"

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
		RgDomain* dom = Domain(i);
		RgSolidDomain* edom = dynamic_cast<RgSolidDomain*>(dom);
		if (edom) edom->BodyForce(R, *this);
	}
}

//-----------------------------------------------------------------------------
// NOTE: Work in progress! Working on integrating body loads as model loads
void FEBodyForce::StiffnessMatrix(FELinearSystem& LS)
{
	for (int i = 0; i<Domains(); ++i)
	{
		RgDomain* dom = Domain(i);
		RgSolidDomain* edom = dynamic_cast<RgSolidDomain*>(dom);
		if (edom) edom->BodyForceStiffness(LS, *this);
	}
}
