#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgPyra13Shape : public RgElementShape
{
public:
    RgPyra13Shape()
        : RgElementShape(ET_PYRA13, 13)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};