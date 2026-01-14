#include "elements/RgElement/RgSolid2dElement.h"

//-----------------------------------------------------------------------------
RgSolid2dElement::RgSolid2dElement(const RgSolid2dElement& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }
}

RgSolid2dElement& RgSolid2dElement::operator=(const RgSolid2dElement& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }
    
    return *this;
}

int RgSolid2dElement::dim()
{
    return 2;
}