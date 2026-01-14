#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>



class RgTet4Shape : public RgSolidElementShape
{
public:
    RgTet4Shape()
        : RgSolidElementShape(ET_TET4, 4)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

