#include "FETet4Shape.h"
#include <vector>
#include <stdexcept>

//=============================================================================
//                                 T E T 4
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
std::vector<double> FETet4::evalH(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet4 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	std::vector<double> H(4);
	H[0] = 1 - r - s - t;
	H[1] = r;
	H[2] = s;
	H[3] = t;
	
	return H;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
std::vector<std::vector<double>> FETet4::evalDeriv(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet4 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return derivatives in the format [dH/dr, dH/ds, dH/dt]
	// Each derivative is a vector of size 4 (one for each node)
	std::vector<std::vector<double>> deriv(3, std::vector<double>(4));

	deriv[0][0] = -1; deriv[1][0] = -1; deriv[2][0] = -1;
	deriv[0][1] = 1;  deriv[1][1] = 0;  deriv[2][1] = 0;
	deriv[0][2] = 0;  deriv[1][2] = 1;  deriv[2][2] = 0;
	deriv[0][3] = 0;  deriv[1][3] = 0;  deriv[2][3] = 1;

	return deriv;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
std::vector<std::vector<double>> FETet4::evalDeriv2(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FETet4 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return second derivatives in the format:
	// [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
	// Each derivative is a vector of size 4 (one for each node)
	std::vector<std::vector<double>> deriv2(6, std::vector<double>(4));

	// All second derivatives are zero for linear shape functions
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			deriv2[i][j] = 0.0;
		}
	}

	return deriv2;
}