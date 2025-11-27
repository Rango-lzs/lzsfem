#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>

class RgHex8Shape : public RgSolidElementShape
{
public:
    RgHex8Shape()
        : RgSolidElementShape(ET_HEX8, 8)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};