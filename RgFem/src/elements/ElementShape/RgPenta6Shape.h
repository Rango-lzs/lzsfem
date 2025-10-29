#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgPenta6Shape : public RgElementShape
{
public:
    RgPenta6Shape()
        : RgElementShape(ET_PENTA6, 6)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};