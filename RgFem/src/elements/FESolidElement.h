#pragma once
#include "elements/RgElement.h"

//-----------------------------------------------------------------------------
//!  This class defines a solid element
//!  对应Abaqus的Continue(Solid)单元，可以是二维的，三维的

class FEM_EXPORT FESolidElement : public FEElement
{
public:
	//! default constructor
	FESolidElement() {}

	//! copy constructor吗
	FESolidElement(const FESolidElement& el);

	//! assignment operator
	FESolidElement& operator = (const FESolidElement& el);

	//! set the element traits
	void SetTraits(FEElementTraits* pt) override;

	double gr(int n) const { return ((FESolidElementTraits*)(m_pTraits))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESolidElementTraits*)(m_pTraits))->gs[n]; }	// integration point coordinate s
	double gt(int n) const { return ((FESolidElementTraits*)(m_pTraits))->gt[n]; }	// integration point coordinate t

	double* GaussWeights() const { return &((FESolidElementTraits*)(m_pTraits))->gw[0]; }			// weights of integration points

	double* Gr(int n) const { return ((FESolidElementTraits*)(m_pTraits))->m_Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESolidElementTraits*)(m_pTraits))->m_Gs[n]; }	// shape function derivative to s
	double* Gt(int n) const { return ((FESolidElementTraits*)(m_pTraits))->m_Gt[n]; }	// shape function derivative to t

	double* Grr(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Grr[n]; }	// shape function 2nd derivative to rr
	double* Gsr(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gsr[n]; }	// shape function 2nd derivative to sr
	double* Gtr(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gtr[n]; }	// shape function 2nd derivative to tr

	double* Grs(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Grs[n]; }	// shape function 2nd derivative to rs
	double* Gss(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gss[n]; }	// shape function 2nd derivative to ss
	double* Gts(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gts[n]; }	// shape function 2nd derivative to ts

	double* Grt(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Grt[n]; }	// shape function 2nd derivative to rt
	double* Gst(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gst[n]; }	// shape function 2nd derivative to st
	double* Gtt(int n) const { return ((FESolidElementTraits*)(m_pTraits))->Gtt[n]; }	// shape function 2nd derivative to tt

																					//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t) const { ((FESolidElementTraits*)(m_pTraits))->shape_fnc(H, r, s, t); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const { ((FESolidElementTraits*)(m_pTraits))->shape_deriv(Hr, Hs, Ht, r, s, t); }

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) const { ((FESolidElementTraits*)(m_pTraits))->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t); }

	Vector3d evaluate(Vector3d* v, double r, double s, double t) const;
	double evaluate(double* v, double r, double s, double t) const;

	double* Gr(int order, int n) const;
	double* Gs(int order, int n) const;
	double* Gt(int order, int n) const;

	void Serialize(DumpStream& ar) override;

public:
	std::vector<bool>    m_bitfc;    //!< flag for interface nodes
	std::vector<Matrix3d>	m_J0i;		//!< inverse of reference Jacobian
};
