
#include "FESolidElement.h"
#include "basicio/DumpStream.h"


//实体单元
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }

	// copy base class data
    m_mat_id = el.m_mat_id;
    m_id = el.m_id;
    m_loc_id = el.m_loc_id;
	m_node = el.m_node;
    m_loc_node = el.m_loc_node;
    /*m_lm = el.m_lm;
    m_val = el.m_val;*/
    m_state = el.m_state;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }

	// copy base class data
    m_mat_id = el.m_mat_id;
	m_id = el.m_id;
	m_loc_id = el.m_loc_id;
	m_node = el.m_node;
	m_loc_node = el.m_loc_node;
    /*m_lm = el.m_lm;
    m_val = el.m_val;*/
	m_state = el.m_state;

	return (*this);
}

void FESolidElement::SetTraits(FEElementTraits* pt)
{
	FEElement::SetTraits(pt);

	int ni = GaussPointSize();
	m_J0i.resize(ni);
}

Vector3d FESolidElement::evaluate(Vector3d* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
    int neln = NodeSize();
	Vector3d p(0, 0, 0);
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double FESolidElement::evaluate(double* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = NodeSize();
	double p = 0.0;
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double* FESolidElement::Gr(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pTraits))->m_Gr_p[order][n] : ((FESolidElementTraits*)(m_pTraits))->m_Gr[n]); }	// shape function derivative to r
double* FESolidElement::Gs(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pTraits))->m_Gs_p[order][n] : ((FESolidElementTraits*)(m_pTraits))->m_Gs[n]); }	// shape function derivative to s
double* FESolidElement::Gt(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pTraits))->m_Gt_p[order][n] : ((FESolidElementTraits*)(m_pTraits))->m_Gt[n]); }	// shape function derivative to t

void FESolidElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		ar & m_J0i;
		ar & m_bitfc;
	}
}
