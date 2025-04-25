#include "elements/FETrussElement.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FETrussElement::FETrussElement()
{
    m_a0 = 0.0;
    m_lam = 1.0;
    m_tau = 0.0;
    m_L0 = 0.0;
}

FETrussElement::FETrussElement(const FETrussElement& el)
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

FETrussElement& FETrussElement::operator=(const FETrussElement& el)
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
void FETrussElement::Serialize(DumpStream& ar)
{
    FEElement::Serialize(ar);
    if (ar.IsShallow() == false)
    {
        ar& m_a0& m_L0& m_lam& m_tau;
    }
}