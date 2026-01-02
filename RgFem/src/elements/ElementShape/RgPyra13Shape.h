#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include <stdexcept>

namespace RgFem {

class RgPyra13Shape : public RgSolidElementShape
{
public:
    RgPyra13Shape()
        : RgElementShape(ET_PYRA13, 13)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem