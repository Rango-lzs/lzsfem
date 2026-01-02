#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"

namespace RgFem {

class RgHex8Shape : public RgSolidElementShape
{
public:
    RgHex8Shape()
        : RgSolidElementShape(ET_HEX8, 8)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem