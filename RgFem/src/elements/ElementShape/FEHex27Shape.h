#pragma once
#include "RgSolidElementShape.h"

//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27 : public RgSolidElementShape
{
public:
	FEHex27() : RgSolidElementShape(ET_HEX27, 27) {}

	std::vector<double> evalH(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;
    std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};