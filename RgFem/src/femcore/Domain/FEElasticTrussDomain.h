#pragma once
#include "femcore/Domain/FETrussDomain.h"
#include "FEElasticDomain.h"
#include "materials/FESolidMaterial.h"
#include "femcore/FEDofList.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FEElasticTrussDomain : public FETrussDomain, public FEElasticDomain
{
public:
	//! Constructor
	FEElasticTrussDomain(FEModel* pfem);

	//! copy operator
	FEElasticTrussDomain& operator = (FEElasticTrussDomain& d);

	//! initialize the domain
	bool Init() override;

	//! Reset data
	void Reset() override;

	//! Initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack truss element data
	void UnpackLM(FEElement& el, std::vector<int>& lm) override;

	//! get the material
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! Activate domain
	void Activate() override;

	//! get the dof list
	const FEDofList& GetDOFList() const override;

public: // overloads from FEElasticDomain

	//! update the truss stresses
	void Update(const FETimeInfo& tp) override;

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! calculate body force \todo implement this
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override { assert(false); }

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override { assert(false); }

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! intertial stiffness matrix
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! body force stiffness matrix \todo implement this
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override { assert(false); }

	//! elemental mass matrix
	void ElementMassMatrix(FETrussElement& el, Matrix& ke);

protected:
	//! calculates the truss element stiffness matrix
	void ElementStiffness(int iel, Matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FETrussElement& el, std::vector<double>& fe);

protected:
	FESolidMaterial*	m_pMat;
	double	m_a0;
	double	m_v;

	FEDofList	m_dofU;

	DECLARE_PARAM_LIST();
};
