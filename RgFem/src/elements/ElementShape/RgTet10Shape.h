#pragma once
#include "elements/ElementShape/RgSolidElementShape.h"
#include <stdexcept>

namespace RgFem {

class RgTet10Shape : public RgSolidElementShape
{
public:
    RgTet10Shape()
        : RgSolidElementShape(ET_TET10, 10)
    {
    }

    std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem