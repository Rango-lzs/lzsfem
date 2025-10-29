#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

class RgTet10Shape : public RgElementShape
{
public:
    RgTet10Shape()
        : RgElementShape(ET_TET10, 10)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};