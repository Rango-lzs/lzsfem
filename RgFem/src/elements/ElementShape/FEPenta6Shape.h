#pragma once
#include "RgSolidElementShape.h"

//=============================================================================
class FEPenta6 : public RgSolidElementShape
{
public:
    FEPenta6()
        : RgSolidElementShape(ET_PENTA6, 6)
    {
    }

    //! values of shape functions
    std::vector<double> evalH(const NaturalCoord& coord) override;

    //! values of shape function derivatives
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;

    //! values of shape function second derivatives
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};