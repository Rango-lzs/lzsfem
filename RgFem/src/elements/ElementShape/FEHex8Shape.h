#pragma once
#include "elements/ElementShape/RgElementShape.h"


class FEHex8 : public RgElementShape
{
public:
    FEHex8()
        : RgElementShape(ET_HEX8, 8)
    {
    }

    std::vector<double> evalH(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv2(NaturalCoord coord) override;
};


//=============================================================================
class FEHex20 : public RgElementShape
{
public:
    FEHex20()
        : RgElementShape(ET_HEX20, 20)
    {
    }

    std::vector<double> evalH(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv2(NaturalCoord coord) override;
};


//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27 : public RgElementShape
{
public:
    FEHex27()
        : RgElementShape(ET_HEX27, 27)
    {
    }

    std::vector<double> evalH(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv(NaturalCoord coord) override;
    std::vector<std::vector<double>> evalDeriv2(NaturalCoord coord) override;
};
