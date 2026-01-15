#pragma once
#include "RgDomain.h"
#include "elements/RgElement/RgTrussElement.h"
#include "femcore/FEDofList.h"

// Forward declaration
class RgAssembler;

//-----------------------------------------------------------------------------
//! Domain described by 3D truss elements
class FEM_EXPORT RgTrussDomain : public RgDomain
{
public:
	//! Constructor
	RgTrussDomain(FEModel* pfem);

	//! copy operator
	RgTrussDomain& operator = (RgTrussDomain& d);

	//! Initialize data
	bool Init() override;

	//! Reset data
	void Reset() override;

	//! Initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

	//! Unpack truss element data
	void UnpackLM(RgElement& el, std::vector<int>& lm) override;

	//! get the material
	RgMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(RgMaterial* pmat) override;

	//! Activate domain
	void Activate() override;

	//! get the dof list
	const FEDofList& GetDOFList() const;

public: // overloads from FEElasticDomain

	//! update the truss stresses
	void Update(const FETimeInfo& tp);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! calculate body force \todo implement this
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) { assert(false); }

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, std::vector<double>& F)  { assert(false); }

	//! calculates the global stiffness Matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) ;

	//! intertial stiffness Matrix
	void MassMatrix(FELinearSystem& LS, double scale) ;

	//! body force stiffness Matrix \todo implement this
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)  { assert(false); }

	//! elemental mass Matrix
	void ElementMassMatrix(RgTrussElement& el, Matrix& ke);

public:
	bool Create(int nsize, FE_Element_Spec espec) override;

	int Elements() const override { return (int)m_Elem.size(); }

	RgTrussElement& Element(int i) { return m_Elem[i]; }

	RgElement& ElementRef(int n) override { return m_Elem[n]; }
	const RgElement& ElementRef(int n) const override { return m_Elem[n]; }

	void ForEachElement(std::function<void(RgElement& el)> f) override;

public:
	void ForEachTrussElement(std::function<void(RgTrussElement& el)> f);

public:
	//! Calculate the truss normal
	Vector3d TrussNormal(RgTrussElement& el);

protected:
	//! calculates the truss element stiffness Matrix
	void ElementStiffness(int iel, Matrix& ke);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(RgTrussElement& el, std::vector<double>& fe);

protected:
	std::vector<RgTrussElement>	m_Elem;
	RgMaterial*	m_pMat;
	double	m_a0;

	FEDofList	m_dofU;
	
	// Assembler for this domain
	RgAssembler* m_assembler;

	DECLARE_PARAM_LIST();
};