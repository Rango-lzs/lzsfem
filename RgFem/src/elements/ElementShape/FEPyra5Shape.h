#pragma once
#include "FESolidElementShape.h"

//=============================================================================
class FEPyra5 : public FESolidElementShape
{
public:
	FEPyra5() : FESolidElementShape(ET_PYRA5, 5) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) override;
};