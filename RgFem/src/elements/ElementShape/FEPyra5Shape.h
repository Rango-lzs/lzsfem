#pragma once
#include "RgSolidElementShape.h"

//=============================================================================
class FEPyra5 : public RgSolidElementShape
{
public:
	FEPyra5() : RgSolidElementShape(ET_PYRA5, 5) {}

	//! values of shape functions
	std::vector<double> evalH(const NaturalCoord& coord) override;

	//! values of shape function derivatives
	std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;

	//! values of shape function second derivatives
	std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override;
};