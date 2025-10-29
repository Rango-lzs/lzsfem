#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgHex27Shape : public RgElementShape
{
public:
    RgHex27Shape()
        : RgElementShape(ET_HEX27, 27)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};