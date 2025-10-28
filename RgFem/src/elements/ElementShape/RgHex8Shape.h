#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgHex8Shape : public RgElementShape
{
public:
    RgHex8Shape()
        : RgElementShape(ET_HEX8, 8)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};