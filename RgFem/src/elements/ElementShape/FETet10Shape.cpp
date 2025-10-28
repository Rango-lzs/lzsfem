#include "FETet10Shape.h"

//=============================================================================
//                                 T E T 1 0
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FETet10::shape_fnc(double* H, double r, double s, double t)
{
	 double u = 1 - r - s - t;
	 
	 H[0] = u*(2*u - 1);
	 H[1] = r*(2*r - 1);
	 H[2] = s*(2*s - 1);
	 H[3] = t*(2*t - 1);
	 H[4] = 4*r*u;
	 H[5] = 4*r*s;
	 H[6] = 4*s*u;
	 H[7] = 4*t*u;
	 H[8] = 4*r*t;
	 H[9] = 4*s*t;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FETet10::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	double u = 1 - r - s - t;

	Hr[0] = 1 - 4 * u;
	Hr[1] = 4 * r - 1;
	Hr[2] = 0;
	Hr[3] = 0;
	Hr[4] = 4 * (u - r);
	Hr[5] = 4 * s;
	Hr[6] = -4 * s;
	Hr[7] = -4 * t;
	Hr[8] = 4 * t;
	Hr[9] = 0;

	Hs[0] = 1 - 4 * u;
	Hs[1] = 0;
	Hs[2] = 4 * s - 1;
	Hs[3] = 0;
	Hs[4] = -4 * r;
	Hs[5] = 4 * r;
	Hs[6] = 4 * (u - s);
	Hs[7] = -4 * t;
	Hs[8] = 0;
	Hs[9] = 4 * t;

	Ht[0] = 1 - 4 * u;
	Ht[1] = 0;
	Ht[2] = 0;
	Ht[3] = 4 * t - 1;
	Ht[4] = -4 * r;
	Ht[5] = 0;
	Ht[6] = -4 * s;
	Ht[7] = 4 * (u - t);
	Ht[8] = 4 * r;
	Ht[9] = 4 * s;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FETet10::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] = 4.0; Hss[0] = 4.0; Htt[0] = 4.0;
	Hrr[1] = 4.0; Hss[1] = 0.0; Htt[1] = 0.0;
	Hrr[2] = 0.0; Hss[2] = 4.0; Htt[2] = 0.0;
	Hrr[3] = 0.0; Hss[3] = 0.0; Htt[3] = 4.0;
	Hrr[4] = -8.0; Hss[4] = 0.0; Htt[4] = 0.0;
	Hrr[5] = 0.0; Hss[5] = 0.0; Htt[5] = 0.0;
	Hrr[6] = 0.0; Hss[6] = -8.0; Htt[6] = 0.0;
	Hrr[7] = 0.0; Hss[7] = 0.0; Htt[7] = -8.0;
	Hrr[8] = 0.0; Hss[8] = 0.0; Htt[8] = 0.0;
	Hrr[9] = 0.0; Hss[9] = 0.0; Htt[9] = 0.0;

	Hrs[0] = 4.0; Hst[0] = 4.0; Hrt[0] = 4.0;
	Hrs[1] = -4.0; Hst[1] = 0.0; Hrt[1] = 0.0;
	Hrs[2] = 0.0; Hst[2] = -4.0; Hrt[2] = 0.0;
	Hrs[3] = 0.0; Hst[3] = 0.0; Hrt[3] = -4.0;
	Hrs[4] = 4.0; Hst[4] = 0.0; Hrt[4] = -4.0;
	Hrs[5] = 4.0; Hst[5] = 4.0; Hrt[5] = 0.0;
	Hrs[6] = -4.0; Hst[6] = 0.0; Hrt[6] = 0.0;
	Hrs[7] = -4.0; Hst[7] = -4.0; Hrt[7] = 0.0;
	Hrs[8] = 0.0; Hst[8] = 0.0; Hrt[8] = 4.0;
	Hrs[9] = 0.0; Hst[9] = 4.0; Hrt[9] = 0.0;
}