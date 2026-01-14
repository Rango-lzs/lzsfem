#include "elements/RgElement/RgTrussElement.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
RgTrussElement::RgTrussElement()
{
    m_a0 = 0.0;
    m_lam = 1.0;
    m_tau = 0.0;
    m_L0 = 0.0;
}

RgTrussElement::RgTrussElement(const RgTrussElement& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }

    // copy base class data
    

    // truss data
    m_a0 = el.m_a0;
    m_L0 = el.m_L0;
    m_lam = el.m_lam;
    m_tau = el.m_tau;
}

RgTrussElement& RgTrussElement::operator=(const RgTrussElement& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }

    // copy base class data
  

    // copy truss data
    m_a0 = el.m_a0;
    m_L0 = el.m_L0;
    m_lam = el.m_lam;
    m_tau = el.m_tau;

    return (*this);
}

//-----------------------------------------------------------------------------
void RgTrussElement::Serialize(DumpStream& ar)
{
    FEElement::Serialize(ar);
    if (ar.IsShallow() == false)
    {
        ar& m_a0& m_L0& m_lam& m_tau;
    }
}