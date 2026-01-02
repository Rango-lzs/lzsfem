#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>

namespace RgFem {

class RgHex27Shape : public RgSolidElementShape
{
public:
    RgHex27Shape()
        : RgSolidElementShape(ET_HEX27, 27)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem