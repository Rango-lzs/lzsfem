#include "FEPenta15Shape.h"
#include <vector>
#include <stdexcept>

//=============================================================================
//                                 P E N T A 1 5
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
std::vector<double> FEPenta15::evalH(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPenta15 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();
	double u = 1 - r - s;

	std::vector<double> H(15);
	H[0] = -((1 - t*t)*u) / 2. + ((1 - t)*u*(-1 + 2 * u)) / 2.;
	H[1] = (r*(-1 + 2 * r)*(1 - t)) / 2. - (r*(1 - t*t)) / 2.;
	H[2] = (s*(-1 + 2 * s)*(1 - t)) / 2. - (s*(1 - t*t)) / 2.;
	H[3] = -((1 - t*t)*u) / 2. + ((1 + t)*u*(-1 + 2 * u)) / 2.;
	H[4] = (r*(-1 + 2 * r)*(1 + t)) / 2. - (r*(1 - t*t)) / 2.;
	H[5] = (s*(-1 + 2 * s)*(1 + t)) / 2. - (s*(1 - t*t)) / 2.;
	H[6] = 2 * r*(1 - t)*u;
	H[7] = 2 * r*s*(1 - t);
	H[8] = 2 * s*(1 - t)*u;
	H[9] = 2 * r*(1 + t)*u;
	H[10] = 2 * r*s*(1 + t);
	H[11] = 2 * s*(1 + t)*u;
	H[12] = (1 - t*t)*u;
	H[13] = r*(1 - t*t);
	H[14] = s*(1 - t*t);
	
	return H;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
std::vector<std::vector<double>> FEPenta15::evalDeriv(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPenta15 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return derivatives in the format [dH/dr, dH/ds, dH/dt]
	// Each derivative is a vector of size 15 (one for each node)
	std::vector<std::vector<double>> deriv(3, std::vector<double>(15));

	deriv[0][0] = -((-1 + t)*(-2 + 4 * r + 4 * s + t)) / 2.;
	deriv[0][1] = (-2 - 4 * r*(-1 + t) + t + t*t) / 2.;
	deriv[0][2] = 0;
	deriv[0][3] = ((-2 + 4 * r + 4 * s - t)*(1 + t)) / 2.;
	deriv[0][4] = ((1 + t)*(-2 + 4 * r + t)) / 2.;
	deriv[0][5] = 0;
	deriv[0][6] = 2 * (-1 + 2 * r + s)*(-1 + t);
	deriv[0][7] = -2 * s*(-1 + t);
	deriv[0][8] = 2 * s*(-1 + t);
	deriv[0][9] = -2 * (-1 + 2 * r + s)*(1 + t);
	deriv[0][10] = 2 * s*(1 + t);
	deriv[0][11] = -2 * s*(1 + t);
	deriv[0][12] = -1 + t*t;
	deriv[0][13] = 1 - t*t;
	deriv[0][14] = 0;

	deriv[1][0] = -((-1 + t)*(-2 + 4 * r + 4 * s + t)) / 2.;
	deriv[1][1] = 0;
	deriv[1][2] = (-2 - 4 * s*(-1 + t) + t + t*t) / 2.;
	deriv[1][3] = ((-2 + 4 * r + 4 * s - t)*(1 + t)) / 2.;
	deriv[1][4] = 0;
	deriv[1][5] = ((1 + t)*(-2 + 4 * s + t)) / 2.;
	deriv[1][6] = 2 * r*(-1 + t);
	deriv[1][7] = -2 * r*(-1 + t);
	deriv[1][8] = 2 * (-1 + r + 2 * s)*(-1 + t);
	deriv[1][9] = -2 * r*(1 + t);
	deriv[1][10] = 2 * r*(1 + t);
	deriv[1][11] = -2 * (-1 + r + 2 * s)*(1 + t);
	deriv[1][12] = -1 + t*t;
	deriv[1][13] = 0;
	deriv[1][14] = 1 - t*t;

	deriv[2][0] = -((-1 + r + s)*(-1 + 2 * r + 2 * s + 2 * t)) / 2.;
	deriv[2][1] = (r*(1 - 2 * r + 2 * t)) / 2.;
	deriv[2][2] = (s*(1 - 2 * s + 2 * t)) / 2.;
	deriv[2][3] = ((-1 + r + s)*(-1 + 2 * r + 2 * s - 2 * t)) / 2.;
	deriv[2][4] = r*(-0.5 + r + t);
	deriv[2][5] = s*(-0.5 + s + t);
	deriv[2][6] = 2 * r*(-1 + r + s);
	deriv[2][7] = -2 * r*s;
	deriv[2][8] = 2 * s*(-1 + r + s);
	deriv[2][9] = -2 * r*(-1 + r + s);
	deriv[2][10] = 2 * r*s;
	deriv[2][11] = -2 * s*(-1 + r + s);
	deriv[2][12] = 2 * (-1 + r + s)*t;
	deriv[2][13] = -2 * r*t;
	deriv[2][14] = -2 * s*t;

	return deriv;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
std::vector<std::vector<double>> FEPenta15::evalDeriv2(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPenta15 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return second derivatives in the format:
	// [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
	// Each derivative is a vector of size 15 (one for each node)
	std::vector<std::vector<double>> deriv2(6, std::vector<double>(15));

	deriv2[0][0] = 2 - 2 * t;
	deriv2[0][1] = 2 - 2 * t;
	deriv2[0][2] = 0;
	deriv2[0][3] = 2 * (1 + t);
	deriv2[0][4] = 2 * (1 + t);
	deriv2[0][5] = 0;
	deriv2[0][6] = 4 * (-1 + t);
	deriv2[0][7] = 0;
	deriv2[0][8] = 0;
	deriv2[0][9] = -4 * (1 + t);
	deriv2[0][10] = 0;
	deriv2[0][11] = 0;
	deriv2[0][12] = 0;
	deriv2[0][13] = 0;
	deriv2[0][14] = 0;

	deriv2[1][0] = 2 - 2 * t;
	deriv2[1][1] = 0;
	deriv2[1][2] = 2 - 2 * t;
	deriv2[1][3] = 2 * (1 + t);
	deriv2[1][4] = 0;
	deriv2[1][5] = 2 * (1 + t);
	deriv2[1][6] = 0;
	deriv2[1][7] = 0;
	deriv2[1][8] = 4 * (-1 + t);
	deriv2[1][9] = 0;
	deriv2[1][10] = 0;
	deriv2[1][11] = -4 * (1 + t);
	deriv2[1][12] = 0;
	deriv2[1][13] = 0;
	deriv2[1][14] = 0;

	deriv2[2][0] = 1 - r - s;
	deriv2[2][1] = r;
	deriv2[2][2] = s;
	deriv2[2][3] = 1 - r - s;
	deriv2[2][4] = r;
	deriv2[2][5] = s;
	deriv2[2][6] = 0;
	deriv2[2][7] = 0;
	deriv2[2][8] = 0;
	deriv2[2][9] = 0;
	deriv2[2][10] = 0;
	deriv2[2][11] = 0;
	deriv2[2][12] = 2 * (-1 + r + s);
	deriv2[2][13] = -2 * r;
	deriv2[2][14] = -2 * s;

	deriv2[3][0] = 2 - 2 * t;
	deriv2[3][1] = 0;
	deriv2[3][2] = 0;
	deriv2[3][3] = 2 * (1 + t);
	deriv2[3][4] = 0;
	deriv2[3][5] = 0;
	deriv2[3][6] = 2 * (-1 + t);
	deriv2[3][7] = 2 - 2 * t;
	deriv2[3][8] = 2 * (-1 + t);
	deriv2[3][9] = -2 * (1 + t);
	deriv2[3][10] = 2 * (1 + t);
	deriv2[3][11] = -2 * (1 + t);
	deriv2[3][12] = 0;
	deriv2[3][13] = 0;
	deriv2[3][14] = 0;

	deriv2[4][0] = 1.5 - 2 * r - 2 * s - t;
	deriv2[4][1] = 0;
	deriv2[4][2] = 0.5 - 2 * s + t;
	deriv2[4][3] = -1.5 + 2 * r + 2 * s - t;
	deriv2[4][4] = 0;
	deriv2[4][5] = -0.5 + 2 * s + t;
	deriv2[4][6] = 2 * r;
	deriv2[4][7] = -2 * r;
	deriv2[4][8] = 2 * (-1 + r + 2 * s);
	deriv2[4][9] = -2 * r;
	deriv2[4][10] = 2 * r;
	deriv2[4][11] = -2 * (-1 + r + 2 * s);
	deriv2[4][12] = 2 * t;
	deriv2[4][13] = 0;
	deriv2[4][14] = -2 * t;

	deriv2[5][0] = 1.5 - 2 * r - 2 * s - t;
	deriv2[5][1] = 0.5 - 2 * r + t;
	deriv2[5][2] = 0;
	deriv2[5][3] = -1.5 + 2 * r + 2 * s - t;
	deriv2[5][4] = -0.5 + 2 * r + t;
	deriv2[5][5] = 0;
	deriv2[5][6] = 2 * (-1 + 2 * r + s);
	deriv2[5][7] = -2 * s;
	deriv2[5][8] = 2 * s;
	deriv2[5][9] = -2 * (-1 + 2 * r + s);
	deriv2[5][10] = 2 * s;
	deriv2[5][11] = -2 * s;
	deriv2[5][12] = 2 * t;
	deriv2[5][13] = -2 * t;
	deriv2[5][14] = 0;

	return deriv2;
}