#pragma once
#include "RgDomain.h"
#include <functional>

//-----------------------------------------------------------------------------
// forward declarations
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;
class FESolidElement;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic domains. 
//! An elastic domain is a domain that can calculate internal forces, stiffness matrices, etc.
class FEM_EXPORT RgElasticDomain
{
public:
	RgElasticDomain(FEModel* pfem);
	virtual ~RgElasticDomain() {}

	//! update stresses
	virtual void Update(const FETimeInfo& tp) = 0;

	//! internal stress forces
	virtual void InternalForces(FEGlobalVector& R) = 0;

	//! calculate body force
	virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf) = 0;

	//! calculates the global stiffness Matrix for this domain
	virtual void StiffnessMatrix(FELinearSystem& LS) = 0;

	//! body force stiffness matrix
	virtual void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) = 0;

	//! mass matrix
	virtual void MassMatrix(FELinearSystem& LS, double scale) = 0;

	//! inertial forces
	virtual void InertialForces(FEGlobalVector& R, std::vector<double>& F) = 0;

	//! loop over all elements
	virtual void ForEachElement(std::function<void(FEElement& el)> f) = 0;

	//! get the material
	virtual FEMaterial* GetMaterial() = 0;

	//! get the dof list
	virtual const FEDofList& GetDOFList() const = 0;

protected:
	FEModel*	m_pfem;
};