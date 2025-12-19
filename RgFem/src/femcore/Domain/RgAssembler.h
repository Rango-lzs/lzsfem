#pragma once
#include "femcore/fem_export.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Abstract interface class for finite element assemblers.
class FEM_EXPORT RgAssembler
{
public:
	RgAssembler(FEModel* pfem);
	virtual ~RgAssembler() {}

	virtual void InternalForces(FEGlobalVector& R) = 0;
	virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf) = 0;
	virtual void InertialForces(FEGlobalVector& R, std::vector<double>& F) = 0;

	//! Calculate global stiffness matrix (only contribution from internal force derivative)
	virtual void StiffnessMatrix(FELinearSystem& LS) = 0;

	//! Calculate stiffness contribution of body forces
	virtual void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) = 0;

	//! Calculate the mass matrix (for dynamic problems)
	virtual void MassMatrix(FELinearSystem& LS, double scale) = 0;

protected:
	FEModel* m_pfem;
};