
#pragma once
#include "femcore/Domain/FEShellDomain.h"
#include "femcore/FEDofList.h"
#include "FEElasticDomain.h"
#include "materials/FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomainOld : public FEShellDomainOld, public FEElasticDomain
{
public:
	FEElasticShellDomainOld(FEModel* pfem);

	//! \todo do I really need this?
	FEElasticShellDomainOld& operator = (FEElasticShellDomainOld& d);

	//! Initialize domain
	bool Init() override;

	//! Activate the domain
	void Activate() override;

	//! Unpack shell element data
	void UnpackLM(FEElement& el, std::vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! get total dof list
	const FEDofList& GetDOFList() const override;

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FESolver* psolver, std::vector<double>& R);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! Calculates inertial forces for dynamic problems | todo implement this (removed assert DSR)
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override { }

	//! calculate body force
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;

	// update stresses
	void Update(const FETimeInfo& tp) override;

	//! calculates the global stiffness Matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	// inertial stiffness \todo implement this (removed assert DSR)
	void MassMatrix(FELinearSystem& LS, double scale) override { }

	// body force stiffness \todo implement this (removed assert DSR)
	void BodyForceStiffness  (FELinearSystem& LS, FEBodyForce& bf) override { }

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElementOld& el, int n, Vector3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElementOld& el, int n, Vector3d g[3]);

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElementOld& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElementOld& el, int n);

    //! calculates covariant basis vectors at an integration point
	void CoBaseVectors(FEShellElementOld& el, int n, Vector3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors(FEShellElementOld& el, int n, Vector3d g[3]);
    
    // jacobian with respect to current frame
	double detJ(FEShellElementOld& el, int n);
    
	// calculate deformation gradient
	double defgrad(FEShellElementOld& el, Matrix3d& F, int n);

	// inverse jacobian with respect to current frame
	double invjact(FEShellElementOld& el, double J[3][3], int n);
 
public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness Matrix
	void ElementStiffness(int iel, Matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress std::vector for shell elements
	void ElementInternalForce(FEShellElementOld& el, std::vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElementOld& el, std::vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElementOld& el, std::vector<double>& fe);

protected:
	FESolidMaterial*	m_pMat;
	FEDofList			m_dofSU;	// shell displacement dofs
	FEDofList			m_dofSR;	// shell rotation dofs
	FEDofList			m_dofR;		// rigid rotation dofs
	FEDofList			m_dof;
};
