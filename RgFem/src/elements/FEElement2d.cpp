#include "elements/FEElement2d.h"

//-----------------------------------------------------------------------------
FEElement2D::FEElement2D(const FEElement2D& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }
}

FEElement2D& FEElement2D::operator=(const FEElement2D& el)
{
    // set the traits of the element
    if (el.m_pTraits)
    {
        SetTraits(el.m_pTraits);
        m_state = el.m_state;
    }

    // copy base class data ??? no need
    /* m_mat = el.m_mat;
     m_nID = el.m_nID;
     m_lid = el.m_lid;
     m_node = el.m_node;
     m_lnode = el.m_lnode;
     m_lm = el.m_lm;
     m_val = el.m_val;*/

    return (*this);
}