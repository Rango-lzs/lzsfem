
//!  This class defines a 2D element
class FEM_EXPORT FEElement2D : public FEElement
{
public:
    //! default constructor
    FEElement2D()
    {
    }

    //! copy constructor
    FEElement2D(const FEElement2D& el);

    //! assignment operator
    FEElement2D& operator=(const FEElement2D& el);

    double* GaussWeights()
    {
        return &((FE2DElementTraits*)(m_pTraits))->gw[0];
    }  // weights of integration points

    double* Hr(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Gr[n];
    }  // shape function derivative to r
    double* Hs(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Gs[n];
    }  // shape function derivative to s

    double* Hrr(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Grr[n];
    }  // shape function 2nd derivative to rr
    double* Hsr(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Gsr[n];
    }  // shape function 2nd derivative to sr

    double* Hrs(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Grs[n];
    }  // shape function 2nd derivative to rs
    double* Hss(int n)
    {
        return ((FE2DElementTraits*)(m_pTraits))->Gss[n];
    }  // shape function 2nd derivative to ss

    //! values of shape functions
    void shape_fnc(double* H, double r, double s)
    {
        ((FE2DElementTraits*)(m_pTraits))->shape(H, r, s);
    }

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s)
    {
        ((FE2DElementTraits*)(m_pTraits))->shape_deriv(Hr, Hs, r, s);
    }
};