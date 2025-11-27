#pragma once
#include "RgSolidElementShape.h"

//=============================================================================
class FETet4 : public RgSolidElementShape
{
public:
	FETet4() : RgSolidElementShape(ET_TET4, 4) {}

	//! values of shape functions
	std::vector<double> evalH(const NaturalCoord& coord) override;

	//! values of shape function derivatives
	std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;

	//! values of shape function second derivatives
	std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};