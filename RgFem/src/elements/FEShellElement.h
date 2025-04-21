#pragma once
#include "RgElement.h"


//-----------------------------------------------------------------------------
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FEM_EXPORT FEShellElement : public FEElement
{
public:
	FEShellElement();

	//! copy constructor
	FEShellElement(const FEShellElement& el);

	//! assignment operator
	FEShellElement& operator = (const FEShellElement& el);

	virtual void SetTraits(FEElementTraits* ptraits) override;

	double gr(int n) { return ((FEShellElementTraits*)(m_pT))->gr[n]; }
	double gs(int n) { return ((FEShellElementTraits*)(m_pT))->gs[n]; }
	double gt(int n) { return ((FEShellElementTraits*)(m_pT))->gt[n]; }

	double* GaussWeights() { return &((FEShellElementTraits*)(m_pT))->gw[0]; }	// weights of integration points

	double* Hr(int n) { return ((FEShellElementTraits*)(m_pT))->Hr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FEShellElementTraits*)(m_pT))->Hs[n]; }	// shape function derivative to s

																			//! values of shape functions
	void shape_fnc(double* H, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_fnc(H, r, s); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_deriv(Hr, Hs, r, s); }

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<double>	m_h0;	//!< initial shell thicknesses
	std::vector<double>	m_ht;	//!< current shell thickness
	std::vector<Vector3d>	m_d0;   //!< initial shell director

	std::vector<Vector3d>	m_g0[3];//!< reference covariant base vectors
	std::vector<Vector3d>	m_gt[3];//!< current covariant base vectors
	std::vector<Vector3d>	m_gp[3];//!< previous covariant base vectors

	std::vector<Vector3d>	m_G0[3];//!< reference contravariant base vectors
	std::vector<Vector3d>	m_Gt[3];//!< current contravariant base vectors

	// indices of solid elements this shell element is attached to.
	// the first element is attached to the back of the shell
	// and the second element is attached to the front.
	// the index is -1 if no solid is attached on that side.
	int        m_elem[2];
};

//-----------------------------------------------------------------------------
// Shell element used by old shell formulation
class FEM_EXPORT FEShellElementOld : public FEShellElement
{
public:
	FEShellElementOld();

	//! copy constructor
	FEShellElementOld(const FEShellElementOld& el);

	//! assignment operator
	FEShellElementOld& operator = (const FEShellElementOld& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<Vector3d>	m_D0;	//!< initial shell directors
};

//-----------------------------------------------------------------------------
// Shell element used by new shell formulations
class FEM_EXPORT FEShellElementNew : public FEShellElement
{
public:
	FEShellElementNew();

	//! copy constructor
	FEShellElementNew(const FEShellElementNew& el);

	//! assignment operator
	FEShellElementNew& operator = (const FEShellElementNew& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public: // EAS parameters

	matrix          m_Kaai;
	matrix          m_fa;
	matrix          m_alpha;
	matrix          m_alphat;
	matrix          m_alphai;
	std::vector<matrix>  m_Kua;
	std::vector<matrix>  m_Kwa;
	std::vector<Matrix3ds>  m_E;
};

