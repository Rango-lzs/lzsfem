#include "FETet4Shape.h"

//=============================================================================
//                                 T E T 4
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FETet4::shape_fnc(double* H, double r, double s, double t)
{
	H[0] = 1 - r - s - t;
	H[1] = r;
	H[2] = s;
	H[3] = t;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FETet4::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[0] = -1; Hs[0] = -1; Ht[0] = -1;
	Hr[1] = 1;	Hs[1] = 0; Ht[1] = 0;
	Hr[2] = 0;	Hs[2] = 1; Ht[2] = 0;
	Hr[3] = 0;	Hs[3] = 0; Ht[3] = 1;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FETet4::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] = 0.0; Hss[0] = 0.0; Htt[0] = 0.0;
	Hrr[1] = 0.0; Hss[1] = 0.0; Htt[1] = 0.0;
	Hrr[2] = 0.0; Hss[2] = 0.0; Htt[2] = 0.0;
	Hrr[3] = 0.0; Hss[3] = 0.0; Htt[3] = 0.0;

	Hrs[0] = 0.0; Hst[0] = 0.0; Hrt[0] = 0.0;
	Hrs[1] = 0.0; Hst[1] = 0.0; Hrt[1] = 0.0;
	Hrs[2] = 0.0; Hst[2] = 0.0; Hrt[2] = 0.0;
	Hrs[3] = 0.0; Hst[3] = 0.0; Hrt[3] = 0.0;
}