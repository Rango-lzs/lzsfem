#include "elements/RgElement/RgSurfaceElement.h"
#include "basicio/DumpStream.h"

//=================================================================================================
// FESurfaceElement
//=================================================================================================

//-----------------------------------------------------------------------------
RgSurfaceElement::RgSurfaceElement()
{
	m_lid = -1;
	m_elem[0] = m_elem[1] = nullptr;
}

RgSurfaceElement::RgSurfaceElement(const RgSurfaceElement& el) : RgElement(el)
{
	// set the traits of the element
	//if (el.m_pTraits) SetTraits(el.m_pTraits);

	// copy base class data

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];
}

RgSurfaceElement& RgSurfaceElement::operator = (const RgSurfaceElement& el)
{
	// make sure the element type is the same
    /*if (m_pTraits == 0) SetTraits(el.m_pTraits);
    else assert(m_pTraits == el.m_pTraits);*/

	// copy base class data

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

