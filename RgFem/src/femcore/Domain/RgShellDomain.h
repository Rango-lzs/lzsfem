#pragma once
#include "RgDomain.h"
#include "elements/FEShellElement.h"

// Forward declaration
class RgAssembler;

//-----------------------------------------------------------------------------
//! Abstract base class for shell domains
class FEM_EXPORT RgShellDomain : public RgDomain
{
    DECLARE_META_CLASS(RgShellDomain, RgDomain);

public:
	RgShellDomain(FEModel* pm);

public:
	virtual bool Create(int nsize, FE_Element_Spec espec) = 0;

	virtual int Elements() const = 0;

	virtual FEElement& ElementRef(int n) = 0;
	virtual const FEElement& ElementRef(int n) const = 0;

	virtual void ForEachElement(std::function<void(FEElement& el)> f);

	virtual void SetMaterial(FEMaterial* pmat);
	
	//! Initialize data
	virtual bool Init();
	
	//! Reset data
	virtual void Reset();
	
	//! Initialize elements
	virtual void PreSolveUpdate(const FETimeInfo& timeInfo);
	
	//! Activate domain
	virtual void Activate();

public:
	//! get the material
	virtual FEMaterial* GetMaterial();

	//! Unpack shell element data
	virtual void UnpackLM(FEElement& el, std::vector<int>& lm);

	//! update the shell stresses
	virtual void Update(const FETimeInfo& tp);

	//! internal stress forces
	virtual void InternalForces(FEGlobalVector& R);

	//! calculate body force
	virtual void BodyForce(FEGlobalVector& R, FEBodyForce& bf);

	//! Calculates inertial forces for dynamic problems
	virtual void InertialForces(FEGlobalVector& R, std::vector<double>& F);

	//! calculates the global stiffness Matrix for this domain
	virtual void StiffnessMatrix(FELinearSystem& LS);

	//! intertial stiffness Matrix
	virtual void MassMatrix(FELinearSystem& LS, double scale);

	//! body force stiffness Matrix
	virtual void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf);

protected:
	// Assembler for this domain
	RgAssembler* m_assembler;
};