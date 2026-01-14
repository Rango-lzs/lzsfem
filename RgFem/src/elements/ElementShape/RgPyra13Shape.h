#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>



class RgPyra13Shape : public RgSolidElementShape
{
public:
    RgPyra13Shape()
        : RgSolidElementShape(ET_PYRA13, 13)
    {
    }

    std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};

