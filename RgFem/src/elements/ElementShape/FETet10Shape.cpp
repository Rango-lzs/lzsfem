#include "FETet10Shape.h"
#include <vector>
#include <stdexcept>

//=============================================================================
//                                 T E T 1 0
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
std::vector<double> FETet10::evalH(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet10 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();
	double u = 1 - r - s - t;

	std::vector<double> H(10);
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
	
	return H;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
std::vector<std::vector<double>> FETet10::evalDeriv(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet10 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();
	double u = 1 - r - s - t;

	// Return derivatives in the format [dH/dr, dH/ds, dH/dt]
	// Each derivative is a vector of size 10 (one for each node)
	std::vector<std::vector<double>> deriv(3, std::vector<double>(10));

	deriv[0][0] = -3.0 + 4.0*r + 4.0*(s + t);
	deriv[0][1] = 4.0*r - 1.0;
	deriv[0][2] = 0.0;
	deriv[0][3] = 0.0;
	deriv[0][4] = 4.0 - 8.0*r - 4.0*(s + t);
	deriv[0][5] = 4.0*s;
	deriv[0][6] = -4.0*s;
	deriv[0][7] = -4.0*t;
	deriv[0][8] = 4.0*t;
	deriv[0][9] = 0.0;

	deriv[1][0] = -3.0 + 4.0*s + 4.0*(r + t);
	deriv[1][1] = 0.0;
	deriv[1][2] = 4.0*s - 1.0;
	deriv[1][3] = 0.0;
	deriv[1][4] = -4.0*r;
	deriv[1][5] = 4.0*r;
	deriv[1][6] = 4.0 - 8.0*s - 4.0*(r + t);
	deriv[1][7] = -4.0*t;
	deriv[1][8] = 0.0;
	deriv[1][9] = 4.0*t;

	deriv[2][0] = -3.0 + 4.0*t + 4.0*(r + s);
	deriv[2][1] = 0.0;
	deriv[2][2] = 0.0;
	deriv[2][3] = 4.0*t - 1.0;
	deriv[2][4] = -4.0*r;
	deriv[2][5] = 0.0;
	deriv[2][6] = -4.0*s;
	deriv[2][7] = 4.0 - 8.0*t - 4.0*(r + s);
	deriv[2][8] = 4.0*r;
	deriv[2][9] = 4.0*s;

	return deriv;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
std::vector<std::vector<double>> FETet10::evalDeriv2(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet10 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return second derivatives in the format:
	// [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
	// Each derivative is a vector of size 10 (one for each node)
	std::vector<std::vector<double>> deriv2(6, std::vector<double>(10));

	deriv2[0][0] = 4.0; deriv2[1][0] = 4.0; deriv2[2][0] = 4.0;
	deriv2[0][1] = 4.0; deriv2[1][1] = 0.0; deriv2[2][1] = 0.0;
	deriv2[0][2] = 0.0; deriv2[1][2] = 4.0; deriv2[2][2] = 0.0;
	deriv2[0][3] = 0.0; deriv2[1][3] = 0.0; deriv2[2][3] = 4.0;
	deriv2[0][4] = -8.0; deriv2[1][4] = 0.0; deriv2[2][4] = 0.0;
	deriv2[0][5] = 0.0; deriv2[1][5] = 0.0; deriv2[2][5] = 0.0;
	deriv2[0][6] = 0.0; deriv2[1][6] = -8.0; deriv2[2][6] = 0.0;
	deriv2[0][7] = 0.0; deriv2[1][7] = 0.0; deriv2[2][7] = -8.0;
	deriv2[0][8] = 0.0; deriv2[1][8] = 0.0; deriv2[2][8] = 0.0;
	deriv2[0][9] = 0.0; deriv2[1][9] = 0.0; deriv2[2][9] = 0.0;

	deriv2[3][0] = 4.0; deriv2[4][0] = 4.0; deriv2[5][0] = 4.0;
	deriv2[3][1] = 0.0; deriv2[4][1] = 0.0; deriv2[5][1] = 0.0;
	deriv2[3][2] = 0.0; deriv2[4][2] = 0.0; deriv2[5][2] = 0.0;
	deriv2[3][3] = 0.0; deriv2[4][3] = 0.0; deriv2[5][3] = 0.0;
	deriv2[3][4] = -4.0; deriv2[4][4] = 0.0; deriv2[5][4] = -4.0;
	deriv2[3][5] = 4.0; deriv2[4][5] = 0.0; deriv2[5][5] = 0.0;
	deriv2[3][6] = -4.0; deriv2[4][6] = -4.0; deriv2[5][6] = 0.0;
	deriv2[3][7] = 0.0; deriv2[4][7] = -4.0; deriv2[5][7] = -4.0;
	deriv2[3][8] = 0.0; deriv2[4][8] = 0.0; deriv2[5][8] = 4.0;
	deriv2[3][9] = 0.0; deriv2[4][9] = 4.0; deriv2[5][9] = 0.0;

	return deriv2;
}