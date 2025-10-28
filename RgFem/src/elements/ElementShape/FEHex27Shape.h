#pragma once
#include "FESolidElementShape.h"

//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27 : public FESolidElementShape
{
public:
	FEHex27() : FESolidElementShape(ET_HEX27, 27) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) override;
};