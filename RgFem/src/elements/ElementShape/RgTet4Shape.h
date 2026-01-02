#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>

namespace RgFem {

class RgTet4Shape : public RgSolidElementShape
{
public:
    RgTet4Shape()
        : RgSolidElementShape(ET_TET4, 4)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem