#pragma once
#include "RgSolidDomain.h"
#include "RgElasticDomain.h"
#include "FEElasticSolidAssembler.h"
#include "materials/FESolidMaterial.h"
#include "femcore/FEDofList.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D solid elements
class FEM_EXPORT RgElasticSolidDomain : public RgSolidDomain, public RgElasticDomain
{
public:
	//! Constructor
	RgElasticSolidDomain(FEModel* pfem);

	//! copy operator
	RgElasticSolidDomain& operator = (RgElasticSolidDomain& d);

	//! Initialize data
	bool Init() override;

	//! Reset data
	void Reset() override;

	//! Initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack solid element data
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

	//! update the solid stresses
	void Update(const FETimeInfo& tp) override;

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! calculate body force \todo implement this
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override { assert(false); }

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, std::vector<double>& F) override { assert(false); }

	//! calculates the global stiffness Matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! intertial stiffness Matrix
	void MassMatrix(FELinearSystem& LS, double scale) override;

	//! body force stiffness Matrix \todo implement this
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override { assert(false); }

	//! elemental mass Matrix
	void ElementMassMatrix(FESolidElement& el, Matrix& ke);

public:
	bool Create(int nsize, FE_Element_Spec espec) override;

	int Elements() const override { return (int)m_Elem.size(); }

	FESolidElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	void ForEachElement(std::function<void(FEElement& el)> f) override;

public:
	void ForEachSolidElement(std::function<void(FESolidElement& el)> f);

protected:
	//! calculates the solid element stiffness Matrix
	void ElementStiffness(int iel, Matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FESolidElement& el, std::vector<double>& fe);

protected:
	std::vector<FESolidElement>	m_Elem;
	FESolidMaterial*	m_pMat;

	FEDofList	m_dofX;
	
	// Assembler for this domain
	FEElasticSolidAssembler* m_assembler;

	DECLARE_PARAM_LIST();
};