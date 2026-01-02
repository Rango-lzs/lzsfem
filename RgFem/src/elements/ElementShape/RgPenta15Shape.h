#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>

namespace RgFem {

class RgPenta15Shape : public RgSolidElementShape
{
public:
    RgPenta15Shape()
        : RgSolidElementShape(ET_PENTA15, 15)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem