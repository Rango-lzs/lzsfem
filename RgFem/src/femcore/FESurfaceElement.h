#pragma once
#include "elements/FEElement.h"

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FEM_EXPORT FESurfaceElement : public FEElement
{
public:
	FESurfaceElement();

	FESurfaceElement(const FESurfaceElement& el);

	FESurfaceElement& operator = (const FESurfaceElement& el);

	virtual void SetTraits(FEElementTraits* pt) override;

	double* GaussWeights() { return &((FESurfaceElementTraits*)(m_pTraits))->gw[0]; }			// weights of integration points
	const double* GaussWeights() const { return &((FESurfaceElementTraits*)(m_pTraits))->gw[0]; }			// weights of integration points
	double gr(int n) const { return ((FESurfaceElementTraits*)(m_pTraits))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESurfaceElementTraits*)(m_pTraits))->gs[n]; }	// integration point coordinate  s
    double cr() const { return ((FESurfaceElementTraits*)(m_pTraits))->cr; }    // centroid point coordinate r
    double cs() const { return ((FESurfaceElementTraits*)(m_pTraits))->cs; }    // centroid point coordinate s

	double* Gr(int n) const { return ((FESurfaceElementTraits*)(m_pTraits))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESurfaceElementTraits*)(m_pTraits))->Gs[n]; }	// shape function derivative to s

	double eval(double* d, int n)
	{
		double* N = H(n);
		int ne = Nodes();
		double a = 0;
		for (int i=0; i<ne; ++i) a += N[i]*d[i];
		return a;
	}

	double eval(int order, double* d, int n)
	{
		double* N = H(order, n);
		int ne = ShapeFunctions(order);
		double a = 0;
		for (int i = 0; i<ne; ++i) a += N[i] * d[i];
		return a;
	}

	double eval(double* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		double a = 0;
		for (int i=0; i<n; ++i) a += H[i]*d[i];
		return a;
	}

	double eval(int order, double* d, double r, double s)
	{
		int n = ShapeFunctions(order);
		double H[FEElement::MAX_NODES];
		shape_fnc(order, H, r, s);
		double a = 0;
		for (int i = 0; i<n; ++i) a += H[i] * d[i];
		return a;
	}

	Vector3d eval(Vector3d* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		Vector3d a(0,0,0);
		for (int i=0; i<n; ++i) a += d[i]*H[i];
		return a;
	}

	Vector3d eval(Vector3d* d, int n)
	{
		int ne = Nodes();
		double* N = H(n);
		Vector3d a(0,0,0);
		for (int i=0; i<ne; ++i) a += d[i]*N[i];
		return a;
	}

	double eval_deriv1(double* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hr[i]*d[i];
		return s;
	}

	double eval_deriv1(int order, double* d, int j)
	{
		double* Hr = Gr(order, j);
		int n = ShapeFunctions(order);
		double s = 0;
		for (int i = 0; i<n; ++i) s += Hr[i] * d[i];
		return s;
	}

	double eval_deriv2(double* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hs[i]*d[i];
		return s;
	}

	Vector3d eval_deriv1(Vector3d* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		Vector3d v(0,0,0);
		for (int i = 0; i<n; ++i) v += d[i]*Hr[i];
		return v;
	}

	Vector3d eval_deriv2(Vector3d* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		Vector3d v(0,0,0);
		for (int i = 0; i<n; ++i) v += d[i]*Hs[i];
		return v;
	}

	double eval_deriv1(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i = 0; i<n; ++i) a += Hr[i] * d[i];
		return a;
	}

	double eval_deriv1(int order, double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(order, Hr, Hs, r, s);
		int n = ShapeFunctions(order);
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hr[i]*d[i];
		return a;
	}

	double eval_deriv2(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hs[i]*d[i];
		return a;
	}

	double eval_deriv2(int order, double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(order, Hr, Hs, r, s);
		int n = ShapeFunctions(order);
		double a = 0;
		for (int i = 0; i<n; ++i) a += Hs[i] * d[i];
		return a;
	}

	void shape_fnc(double* H, double r, double s)
	{
		((FESurfaceElementTraits*)m_pTraits)->shape_fnc(H, r, s);
	}

	void shape_deriv(double* Gr, double* Gs, double r, double s)
	{
		((FESurfaceElementTraits*)m_pTraits)->shape_deriv(Gr, Gs, r, s);
	}

	void shape_deriv(int order, double* Gr, double* Gs, double r, double s)
	{
		((FESurfaceElementTraits*)m_pTraits)->shape_deriv(order, Gr, Gs, r, s);
	}

	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
	{
		((FESurfaceElementTraits*)m_pTraits)->shape_deriv2(Grr, Grs, Gss, r, s);
	}

	void shape_fnc(int order, double* H, double r, double s)
	{
		((FESurfaceElementTraits*)m_pTraits)->shape_fnc(order, H, r, s);
	}

    //! return number of edges
    int facet_edges() const;
    
    //! return node list of edge
    void facet_edge(int j, int* en) const;

	//! serialize
	void Serialize(DumpStream& ar) override;

	double* Gr(int order, int n) const;
	double* Gs(int order, int n) const;
    
public:
    //! local ID of surface element
	int		m_lid;

	// indices of solid or shell element this surface is a face of
	// For solids, a surface element can be connected to two elements 
	// if the surface is an inside surface. For boundary surfaces
	// the second element index is -1. 
	FEElement*		m_elem[2];
};

