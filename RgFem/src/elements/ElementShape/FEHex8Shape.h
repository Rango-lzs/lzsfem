#pragma once
#include "FESolidElementShape.h"

//=============================================================================
class FEHex8 : public FESolidElementShape
{
public:
	FEHex8() : FESolidElementShape(ET_HEX8, 8) {}

	//! values of shape functions
	//! H<8>
	void shape_fnc(double* H, double r, double s, double t) override;

	//! values of shape function derivatives
	//! H<3,8>
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) override;

	//! values of shape function second derivatives
	//! H<6,8>
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) override;
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
