#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgPenta15Shape : public RgElementShape
{
public:
    RgPenta15Shape()
        : RgElementShape(ET_PENTA15, 15)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};