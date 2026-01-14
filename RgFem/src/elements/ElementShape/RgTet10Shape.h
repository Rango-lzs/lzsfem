#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>



class RgTet10Shape : public RgSolidElementShape
{
public:
    RgTet10Shape()
        : RgSolidElementShape(ET_TET10, 10)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

