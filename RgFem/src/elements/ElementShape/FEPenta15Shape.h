#pragma once
#include "FESolidElementShape.h"

//=============================================================================
class FEPenta15 : public FESolidElementShape
{
public:
    FEPenta15()
        : FESolidElementShape(ET_PENTA15, 15)
    {
    }

    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t) override;

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) override;

    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t) override;
};