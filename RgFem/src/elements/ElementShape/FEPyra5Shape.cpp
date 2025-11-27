#include "FEPyra5Shape.h"
#include <vector>
#include <stdexcept>

//=============================================================================
//              P Y R A 5
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
std::vector<double> FEPyra5::evalH(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPyra5 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	std::vector<double> H(5);
	H[0] = 0.125*(1.0 - r)*(1.0 - s)*(1.0 - t);
	H[1] = 0.125*(1.0 + r)*(1.0 - s)*(1.0 - t);
	H[2] = 0.125*(1.0 + r)*(1.0 + s)*(1.0 - t);
	H[3] = 0.125*(1.0 - r)*(1.0 + s)*(1.0 - t);
	H[4] = 0.5*(1.0 + t);
	
	return H;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
std::vector<std::vector<double>> FEPyra5::evalDeriv(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPyra5 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return derivatives in the format [dH/dr, dH/ds, dH/dt]
	// Each derivative is a vector of size 5 (one for each node)
	std::vector<std::vector<double>> deriv(3, std::vector<double>(5));

	deriv[0][0] = -0.125*(1.0 - s)*(1.0 - t);
	deriv[0][1] = 0.125*(1.0 - s)*(1.0 - t);
	deriv[0][2] = 0.125*(1.0 + s)*(1.0 - t);
	deriv[0][3] = -0.125*(1.0 + s)*(1.0 - t);
	deriv[0][4] = 0.0;

	deriv[1][0] = -0.125*(1.0 - r)*(1.0 - t);
	deriv[1][1] = -0.125*(1.0 + r)*(1.0 - t);
	deriv[1][2] = 0.125*(1.0 + r)*(1.0 - t);
	deriv[1][3] = 0.125*(1.0 - r)*(1.0 - t);
	deriv[1][4] = 0.0;

	deriv[2][0] = -0.125*(1.0 - r)*(1.0 - s);
	deriv[2][1] = -0.125*(1.0 + r)*(1.0 - s);
	deriv[2][2] = -0.125*(1.0 + r)*(1.0 + s);
	deriv[2][3] = -0.125*(1.0 - r)*(1.0 + s);
	deriv[2][4] = 0.5;

	return deriv;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
std::vector<std::vector<double>> FEPyra5::evalDeriv2(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEPyra5 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return second derivatives in the format:
	// [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
	// Each derivative is a vector of size 5 (one for each node)
	std::vector<std::vector<double>> deriv2(6, std::vector<double>(5));

	for (int i = 0; i < 5; i++) {
		deriv2[0][i] = 0.0;
		deriv2[1][i] = 0.0;
		deriv2[2][i] = 0.0;
	}

	deriv2[3][0] = 0.125*(1.0 - t); deriv2[4][0] = 0.125*(1.0 - r); deriv2[5][0] = 0.125*(1.0 - s);
	deriv2[3][1] = -0.125*(1.0 - t); deriv2[4][1] = 0.125*(1.0 + r); deriv2[5][1] = -0.125*(1.0 - s);
	deriv2[3][2] = 0.125*(1.0 - t); deriv2[4][2] = -0.125*(1.0 + r); deriv2[5][2] = -0.125*(1.0 + s);
	deriv2[3][3] = -0.125*(1.0 - t); deriv2[4][3] = -0.125*(1.0 - r); deriv2[5][3] = 0.125*(1.0 + s);
	deriv2[3][4] = 0.0; deriv2[4][4] = 0.0; deriv2[5][4] = 0.0;

	return deriv2;
}