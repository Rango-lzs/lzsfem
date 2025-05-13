#pragma once
#include "femcore/Domain/FEDomain.h"
#include "elements/FEShellElement.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell element domains
class FEM_EXPORT FEShellDomain : public FEDomain
{
	DECLARE_META_CLASS(FEShellDomain, FEDomain);

public:
	//! constructor
	FEShellDomain(FEModel* fem);

	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! Reset element data
	void Reset();

	// get a shell element
	virtual FEShellElement& Element(int i) = 0;

	// get the element type (TODO: Move to FEDomain class?)
	int GetElementType() { /*return ElementRef(0).Type();*/ return 0; };

public:
	// evaluate volume of element in reference frame
	virtual double Volume(FEShellElement& el) { return 0.0; }

	// evaluate volume of element in current frame
	virtual double CurrentVolume(FEShellElement& el) { return 0.0; }

	// Initialize shell data (Called from FEMesh::InitShells)
	virtual void InitShells();

public:
    //! get the current nodal coordinates
    void GetCurrentNodalCoordinates(const FEShellElement& el, Vector3d* rt, const bool back = false);
    void GetCurrentNodalCoordinates(const FEShellElement& el, Vector3d* rt, double alpha, const bool back = false);
    
    //! get the reference nodal coordinates
    void GetReferenceNodalCoordinates(const FEShellElement& el, Vector3d* r0, const bool back = false);
    
    //! get the nodal coordinates at previous state
    void GetPreviousNodalCoordinates(const FEShellElement& el, Vector3d* rp, const bool back = false);
    
public:
	void ForEachShellElement(std::function<void(FEShellElement& el)> f);
};

//-----------------------------------------------------------------------------
// Old director-based shell formulation
class FEM_EXPORT FEShellDomainOld : public FEShellDomain
{
public:
	FEShellDomainOld(FEModel* fem);

	//! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEShellElementOld& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

	void InitShells() override;

protected:
	std::vector<FEShellElementOld>	m_Elem;	//!< array of elements
};

//-----------------------------------------------------------------------------
// New shell formulation
class FEM_EXPORT FEShellDomainNew : public FEShellDomain
{
public:
	FEShellDomainNew(FEModel* fem);

	//! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

public:
	//! return nr of elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) override { return m_Elem[n]; }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEShellElementNew& ShellElement(int i) { return m_Elem[i]; }

	double Volume(FEShellElement& el) override;

	double DefaultShellThickness() const { return m_h0; }

	void AssignDefaultShellThickness();

protected:
	double	m_h0;

protected:
	std::vector<FEShellElementNew>	m_Elem;	//!< array of elements

	DECLARE_PARAM_LIST();
};
