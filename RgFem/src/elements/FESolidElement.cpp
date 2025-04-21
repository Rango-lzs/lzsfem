
#include "FESolidElement.h"
#include "DumpStream.h"


//实体单元
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
	m_status = el.m_status;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
	m_status = el.m_status;

	return (*this);
}

void FESolidElement::SetTraits(FEElementTraits* pt)
{
	FEElement::SetTraits(pt);

	int ni = GaussPoints();
	m_J0i.resize(ni);
}

Vector3d FESolidElement::evaluate(Vector3d* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = Nodes();
	Vector3d p(0, 0, 0);
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double FESolidElement::evaluate(double* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = Nodes();
	double p = 0.0;
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double* FESolidElement::Gr(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gr_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gr[n]); }	// shape function derivative to r
double* FESolidElement::Gs(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gs_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gs[n]); }	// shape function derivative to s
double* FESolidElement::Gt(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gt_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gt[n]); }	// shape function derivative to t

void FESolidElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		ar & m_J0i;
		ar & m_bitfc;
	}
}
