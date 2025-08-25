#include "elements/FESurfaceElement.h"
#include "basicio/DumpStream.h"

//=================================================================================================
// FESurfaceElement
//=================================================================================================

//-----------------------------------------------------------------------------
FESurfaceElement::FESurfaceElement()
{
	m_lid = -1;
	m_elem[0] = m_elem[1] = nullptr;
}

FESurfaceElement::FESurfaceElement(const FESurfaceElement& el) : FEElement(el)
{
	// set the traits of the element
	if (el.m_pTraits) SetTraits(el.m_pTraits);

	// copy base class data

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];
}

FESurfaceElement& FESurfaceElement::operator = (const FESurfaceElement& el)
{
	// make sure the element type is the same
	if (m_pTraits == 0) SetTraits(el.m_pTraits);
	else assert(m_pTraits == el.m_pTraits);

	// copy base class data

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

void FESurfaceElement::SetTraits(FEElementTraits* pt)
{
	m_pTraits = pt;
	m_node.resize(NodeSize());
	m_loc_node.resize(NodeSize());
	m_state.Create(GaussPointSize());
}

int FESurfaceElement::facet_edges() const
{
	int nn = NodeSize(), nf = 0;
	switch (nn)
	{
	case 3:
	case 6:
	case 7:
		nf = 3;
		break;
	case 4:
	case 8:
	case 9:
		nf = 4;
		break;
	default:
		assert(false);
	}
	return nf;
}

void FESurfaceElement::facet_edge(int j, int* en) const
{
	int nn = NodeSize();
	switch (nn)
	{
	case 3:
		en[0] = m_loc_node[j];
		en[1] = m_loc_node[(j + 1) % 3];
		break;
	case 6:
	case 7:
		en[0] = m_loc_node[j];
		en[1] = m_loc_node[j + 3];
		en[2] = m_loc_node[(j + 1) % 3];
		break;
	case 4:
		en[0] = m_loc_node[j];
		en[1] = m_loc_node[(j + 1) % 4];
		break;
	case 8:
	case 9:
		en[0] = m_loc_node[j];
		en[1] = m_loc_node[j + 4];
		en[2] = m_loc_node[(j + 1) % 4];
		break;
	}
}

double* FESurfaceElement::Gr(int order, int n) const { return (order >= 0 ? ((FESurfaceElementTraits*)(m_pTraits))->Gr_p[order][n] : ((FESurfaceElementTraits*)(m_pTraits))->Gr[n]); }	// shape function derivative to r
double* FESurfaceElement::Gs(int order, int n) const { return (order >= 0 ? ((FESurfaceElementTraits*)(m_pTraits))->Gs_p[order][n] : ((FESurfaceElementTraits*)(m_pTraits))->Gs[n]); }	// shape function derivative to s

void FESurfaceElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_lid;
	// TODO: Serialize m_elem
}
