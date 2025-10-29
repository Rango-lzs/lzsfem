#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgPyra5Shape : public RgElementShape
{
public:
    RgPyra5Shape()
        : RgElementShape(ET_PYRA5, 5)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};