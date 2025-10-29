#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgTet4Shape : public RgElementShape
{
public:
    RgTet4Shape()
        : RgElementShape(ET_TET4, 4)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};