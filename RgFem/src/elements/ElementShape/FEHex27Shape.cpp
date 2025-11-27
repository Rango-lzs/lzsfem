#include "FEHex27Shape.h"
#include <vector>
#include <stdexcept>

//=============================================================================
//              H E X 2 7
//=============================================================================
// Lookup table for 27-node hex that maps a node index into a triplet that identifies
// the shape functions.
static int HEX27_LUT[27][3] = {
	{ 0,0,0 },
	{ 1,0,0 },
	{ 1,1,0 },
	{ 0,1,0 },
	{ 0,0,1 },
	{ 1,0,1 },
	{ 1,1,1 },
	{ 0,1,1 },
	{ 2,0,0 },
	{ 1,2,0 },
	{ 2,1,0 },
	{ 0,2,0 },
	{ 2,0,1 },
	{ 1,2,1 },
	{ 2,1,1 },
	{ 0,2,1 },
	{ 0,0,2 },
	{ 1,0,2 },
	{ 1,1,2 },
	{ 0,1,2 },
	{ 2,0,2 },
	{ 1,2,2 },
	{ 2,1,2 },
	{ 0,2,2 },
	{ 2,2,0 },
	{ 2,2,1 },
	{ 2,2,2 }
};

//-----------------------------------------------------------------------------
//! values of shape functions
std::vector<double> FEHex27::evalH(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEHex27 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	std::vector<double> H(27);

	double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
	double T[3] = { 0.5*t*(t - 1.0), 0.5*t*(t + 1.0), 1.0 - t*t };

	H[0] = R[0] * S[0] * T[0];
	H[1] = R[1] * S[0] * T[0];
	H[2] = R[1] * S[1] * T[0];
	H[3] = R[0] * S[1] * T[0];
	H[4] = R[0] * S[0] * T[1];
	H[5] = R[1] * S[0] * T[1];
	H[6] = R[1] * S[1] * T[1];
	H[7] = R[0] * S[1] * T[1];
	H[8] = R[2] * S[0] * T[0];
	H[9] = R[1] * S[2] * T[0];
	H[10] = R[2] * S[1] * T[0];
	H[11] = R[0] * S[2] * T[0];
	H[12] = R[2] * S[0] * T[1];
	H[13] = R[1] * S[2] * T[1];
	H[14] = R[2] * S[1] * T[1];
	H[15] = R[0] * S[2] * T[1];
	H[16] = R[0] * S[0] * T[2];
	H[17] = R[1] * S[0] * T[2];
	H[18] = R[1] * S[1] * T[2];
	H[19] = R[0] * S[1] * T[2];
	H[20] = R[2] * S[0] * T[2];
	H[21] = R[1] * S[2] * T[2];
	H[22] = R[2] * S[1] * T[2];
	H[23] = R[0] * S[2] * T[2];
	H[24] = R[2] * S[2] * T[0];
	H[25] = R[2] * S[2] * T[1];
	H[26] = R[2] * S[2] * T[2];

	return H;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
std::vector<std::vector<double>> FEHex27::evalDeriv(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEHex27 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return derivatives in the format [dH/dr, dH/ds, dH/dt]
	// Each derivative is a vector of size 27 (one for each node)
	std::vector<std::vector<double>> deriv(3, std::vector<double>(27));

	double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
	double T[3] = { 0.5*t*(t - 1.0), 0.5*t*(t + 1.0), 1.0 - t*t };

	double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
	double DS[3] = { s - 0.5, s + 0.5, -2.0*s };
	double DT[3] = { t - 0.5, t + 0.5, -2.0*t };

	deriv[0][0] = DR[0] * S[0] * T[0]; deriv[1][0] = R[0] * DS[0] * T[0]; deriv[2][0] = R[0] * S[0] * DT[0];
	deriv[0][1] = DR[1] * S[0] * T[0]; deriv[1][1] = R[1] * DS[0] * T[0]; deriv[2][1] = R[1] * S[0] * DT[0];
	deriv[0][2] = DR[1] * S[1] * T[0]; deriv[1][2] = R[1] * DS[1] * T[0]; deriv[2][2] = R[1] * S[1] * DT[0];
	deriv[0][3] = DR[0] * S[1] * T[0]; deriv[1][3] = R[0] * DS[1] * T[0]; deriv[2][3] = R[0] * S[1] * DT[0];
	deriv[0][4] = DR[0] * S[0] * T[1]; deriv[1][4] = R[0] * DS[0] * T[1]; deriv[2][4] = R[0] * S[0] * DT[1];
	deriv[0][5] = DR[1] * S[0] * T[1]; deriv[1][5] = R[1] * DS[0] * T[1]; deriv[2][5] = R[1] * S[0] * DT[1];
	deriv[0][6] = DR[1] * S[1] * T[1]; deriv[1][6] = R[1] * DS[1] * T[1]; deriv[2][6] = R[1] * S[1] * DT[1];
	deriv[0][7] = DR[0] * S[1] * T[1]; deriv[1][7] = R[0] * DS[1] * T[1]; deriv[2][7] = R[0] * S[1] * DT[1];
	deriv[0][8] = DR[2] * S[0] * T[0]; deriv[1][8] = R[2] * DS[0] * T[0]; deriv[2][8] = R[2] * S[0] * DT[0];
	deriv[0][9] = DR[1] * S[2] * T[0]; deriv[1][9] = R[1] * DS[2] * T[0]; deriv[2][9] = R[1] * S[2] * DT[0];
	deriv[0][10] = DR[2] * S[1] * T[0]; deriv[1][10] = R[2] * DS[1] * T[0]; deriv[2][10] = R[2] * S[1] * DT[0];
	deriv[0][11] = DR[0] * S[2] * T[0]; deriv[1][11] = R[0] * DS[2] * T[0]; deriv[2][11] = R[0] * S[2] * DT[0];
	deriv[0][12] = DR[2] * S[0] * T[1]; deriv[1][12] = R[2] * DS[0] * T[1]; deriv[2][12] = R[2] * S[0] * DT[1];
	deriv[0][13] = DR[1] * S[2] * T[1]; deriv[1][13] = R[1] * DS[2] * T[1]; deriv[2][13] = R[1] * S[2] * DT[1];
	deriv[0][14] = DR[2] * S[1] * T[1]; deriv[1][14] = R[2] * DS[1] * T[1]; deriv[2][14] = R[2] * S[1] * DT[1];
	deriv[0][15] = DR[0] * S[2] * T[1]; deriv[1][15] = R[0] * DS[2] * T[1]; deriv[2][15] = R[0] * S[2] * DT[1];
	deriv[0][16] = DR[0] * S[0] * T[2]; deriv[1][16] = R[0] * DS[0] * T[2]; deriv[2][16] = R[0] * S[0] * DT[2];
	deriv[0][17] = DR[1] * S[0] * T[2]; deriv[1][17] = R[1] * DS[0] * T[2]; deriv[2][17] = R[1] * S[0] * DT[2];
	deriv[0][18] = DR[1] * S[1] * T[2]; deriv[1][18] = R[1] * DS[1] * T[2]; deriv[2][18] = R[1] * S[1] * DT[2];
	deriv[0][19] = DR[0] * S[1] * T[2]; deriv[1][19] = R[0] * DS[1] * T[2]; deriv[2][19] = R[0] * S[1] * DT[2];
	deriv[0][20] = DR[2] * S[0] * T[2]; deriv[1][20] = R[2] * DS[0] * T[2]; deriv[2][20] = R[2] * S[0] * DT[2];
	deriv[0][21] = DR[1] * S[2] * T[2]; deriv[1][21] = R[1] * DS[2] * T[2]; deriv[2][21] = R[1] * S[2] * DT[2];
	deriv[0][22] = DR[2] * S[1] * T[2]; deriv[1][22] = R[2] * DS[1] * T[2]; deriv[2][22] = R[2] * S[1] * DT[2];
	deriv[0][23] = DR[0] * S[2] * T[2]; deriv[1][23] = R[0] * DS[2] * T[2]; deriv[2][23] = R[0] * S[2] * DT[2];
	deriv[0][24] = DR[2] * S[2] * T[0]; deriv[1][24] = R[2] * DS[2] * T[0]; deriv[2][24] = R[2] * S[2] * DT[0];
	deriv[0][25] = DR[2] * S[2] * T[1]; deriv[1][25] = R[2] * DS[2] * T[1]; deriv[2][25] = R[2] * S[2] * DT[1];
	deriv[0][26] = DR[2] * S[2] * T[2]; deriv[1][26] = R[2] * DS[2] * T[2]; deriv[2][26] = R[2] * S[2] * DT[2];

	return deriv;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
std::vector<std::vector<double>> FEHex27::evalDeriv2(const NaturalCoord& coord)
{
	const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
	if (!coord3d) {
		throw std::invalid_argument("FEHex27 requires NaturalCoord3d coordinates");
	}

	double r = coord3d->getR();
	double s = coord3d->getS();
	double t = coord3d->getT();

	// Return second derivatives in the format:
	// [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
	// Each derivative is a vector of size 27 (one for each node)
	std::vector<std::vector<double>> deriv2(6, std::vector<double>(27));

	double NR[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double NS[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
	double NT[3] = { 0.5*t*(t - 1.0), 0.5*t*(t + 1.0), 1.0 - t*t };

	double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
	double DS[3] = { s - 0.5, s + 0.5, -2.0*s };
	double DT[3] = { t - 0.5, t + 0.5, -2.0*t };

	double HR[3] = { 1.0, 1.0, -2.0 };
	double HS[3] = { 1.0, 1.0, -2.0 };
	double HT[3] = { 1.0, 1.0, -2.0 };

	for (int a = 0; a < 27; ++a)
	{
		int i = HEX27_LUT[a][0];
		int j = HEX27_LUT[a][1];
		int k = HEX27_LUT[a][2];

		deriv2[0][a] = HR[i] * NS[j] * NT[k];  // Hrr
		deriv2[1][a] = NR[i] * HS[j] * NT[k];  // Hss
		deriv2[2][a] = NR[i] * NS[j] * HT[k];  // Htt
		deriv2[3][a] = DR[i] * DS[j] * NT[k];  // Hrs
		deriv2[4][a] = NR[i] * DS[j] * DT[k];  // Hst
		deriv2[5][a] = DR[i] * NS[j] * DT[k];  // Hrt
	}

	return deriv2;
}