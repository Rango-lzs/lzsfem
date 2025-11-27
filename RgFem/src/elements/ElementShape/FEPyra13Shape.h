#pragma once
#include "RgSolidElementShape.h"

//=============================================================================
class FEPyra13 : public RgSolidElementShape
{
public:
    FEPyra13() : RgSolidElementShape(ET_PYRA13, 13) {}
    
    //! values of shape functions
    std::vector<double> evalH(const NaturalCoord& coord) override;
    
    //! values of shape function derivatives
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    
    //! values of shape function second derivatives
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};