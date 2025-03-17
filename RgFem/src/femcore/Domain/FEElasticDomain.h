#pragma once
#include "femcore/fem_export.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FESolver;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic domains.
// 1¡¢×é×°¾ØÕó
// 2¡¢¼ÆËã²Ð²î
class FEM_EXPORT FEElasticDomain
{
public:
	FEElasticDomain(FEModel* pfem);
	virtual ~FEElasticDomain(){}


	virtual void InternalForces(FEGlobalVector& R) = 0;
	virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf) = 0;
	virtual void InertialForces(FEGlobalVector& R, std::vector<double>& F) = 0;



	//! Calculate global stiffness matrix (only contribution from internal force derivative)
	//! \todo maybe I should rename this the InternalStiffness matrix?
	virtual void StiffnessMatrix   (FELinearSystem& LS) = 0;

	//! Calculate stiffness contribution of body forces
	virtual void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) = 0;

	//! calculate the mass matrix (for dynamic problems)
	virtual void MassMatrix(FELinearSystem& LS, double scale) = 0;
};
