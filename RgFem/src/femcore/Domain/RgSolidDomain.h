#pragma once
#include "RgDomain.h"
#include "elements/RgElement/RgSolidElement.h"

// Forward declaration
class RgAssembler;

//-----------------------------------------------------------------------------
//! Abstract base class for solid domains
class FEM_EXPORT RgSolidDomain : public RgDomain
{
    DECLARE_META_CLASS(RgSolidDomain, RgDomain);

public:
	RgSolidDomain();
	RgSolidDomain(FEModel* pm);
	~RgSolidDomain();

public:
	virtual bool Create(int nsize, ElementType espec);
	RgSolidElement* CreateSolidElement(ElementType type);

	virtual int Elements() const;

	virtual RgElement& ElementRef(int n);
	virtual const RgElement& ElementRef(int n) const;

	virtual void ForEachElement(std::function<void(RgElement& el)> f);

	virtual void SetMaterial(RgMaterial* pmat);
	
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
	virtual RgMaterial* GetMaterial();

	//! Unpack solid element data
	virtual void UnpackLM(RgElement& el, std::vector<int>& lm);

	//! update the solid stresses
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