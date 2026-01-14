#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>


class RgHex20Shape : public RgSolidElementShape
{
public:
    RgHex20Shape()
        : RgSolidElementShape(ET_HEX20, 20)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

