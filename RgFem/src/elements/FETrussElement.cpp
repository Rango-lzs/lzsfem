//=================================================================================================
// FETrussElement
//=================================================================================================

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
    m_mat = el.m_mat;
    m_nID = el.m_nID;
    m_lid = el.m_lid;
    m_node = el.m_node;
    m_lnode = el.m_lnode;
    m_lm = el.m_lm;
    m_val = el.m_val;

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
    m_mat = el.m_mat;
    m_nID = el.m_nID;
    m_lid = el.m_lid;
    m_node = el.m_node;
    m_lnode = el.m_lnode;
    m_lm = el.m_lm;
    m_val = el.m_val;

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