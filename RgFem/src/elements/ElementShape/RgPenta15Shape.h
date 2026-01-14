#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>



class RgPenta15Shape : public RgSolidElementShape
{
public:
    RgPenta15Shape()
        : RgSolidElementShape(ET_PENTA15, 15)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

