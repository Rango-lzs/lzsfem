#pragma once
#include "elements/RgElement/RgElement.h"
#include <vector>

//! This class defines a solid element
class RgSolidElementTraits;

//! Corresponds to Abaqus Continuum(Solid) elements, which are multidimensional (3D)
class FEM_EXPORT RgSolidElement : public RgElement
{
public:
	//! default constructor
	RgSolidElement() = default;
	~RgSolidElement() = default;

	//! copy constructor
	RgSolidElement(const RgSolidElement& el);

	//! assignment operator
	RgSolidElement& operator = (const RgSolidElement& el);

	//! set the element traits
	void SetTraits(FEElementTraits* pt) override;

	//n : the n-th gauss point
	double gr(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->gs[n]; }	// integration point coordinate s
	double gt(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->gt[n]; }	// integration point coordinate t

	// returns weights of integration points as a vector
	std::vector<double> GaussWeights() const { return ((RgSolidElementTraits*)(m_pTraits))->gw; }

	// dH/dr[n] -- return shape function derivative arrays as vectors
	std::vector<double> Gr(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->m_Gr[n]; }	// shape function derivative to r
	std::vector<double> Gs(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->m_Gs[n]; }	// shape function derivative to s
	std::vector<double> Gt(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->m_Gt[n]; }	// shape function derivative to t

	// dH2/dr2[n] -- second derivatives as vectors
	std::vector<double> Grr(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Grr[n]; }	// shape function 2nd derivative to rr
	std::vector<double> Gsr(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gsr[n]; }	// shape function 2nd derivative to sr
	std::vector<double> Gtr(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gtr[n]; }	// shape function 2nd derivative to tr

	std::vector<double> Grs(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Grs[n]; }	// shape function 2nd derivative to rs
	std::vector<double> Gss(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gss[n]; }	// shape function 2nd derivative to ss
	std::vector<double> Gts(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gts[n]; }	// shape function 2nd derivative to ts

	std::vector<double> Grt(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Grt[n]; }	// shape function 2nd derivative to rt
	std::vector<double> Gst(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gst[n]; }	// shape function 2nd derivative to st
	std::vector<double> Gtt(int n) const { return ((RgSolidElementTraits*)(m_pTraits))->Gtt[n]; }	// shape function 2nd derivative to tt

	//! values of shape functions (unchanged API)
	void shape_fnc(double* H, double r, double s, double t) const { ((RgSolidElementTraits*)(m_pTraits))->shape_fnc(H, r, s, t); }

	//! values of shape function derivatives (unchanged API)
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const { ((RgSolidElementTraits*)(m_pTraits))->shape_deriv(Hr, Hs, Ht, r, s, t); }

	//! values of shape function second derivatives (unchanged API)
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) const { ((RgSolidElementTraits*)(m_pTraits))->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t); }

	Vector3d evaluate(Vector3d* v, double r, double s, double t) const;
	double evaluate(double* v, double r, double s, double t) const;

	// overloads that select derivative order (now return vectors)
	std::vector<double> Gr(int order, int n) const;
	std::vector<double> Gs(int order, int n) const;
	std::vector<double> Gt(int order, int n) const;

	void Serialize(DumpStream& ar) override;

	virtual void StiffnessMatrix(Matrix& ke);
	virtual void MaterialStiffness(Matrix& ke);
	virtual void GeometricalStiffness(Matrix& ke);


public:
	std::vector<bool>    m_bitfc;    //!< flag for interface nodes
	std::vector<Matrix3d>	m_J0i;		//!< inverse of reference Jacobian
};
