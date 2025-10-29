#include "elements/ElementShape/RgPenta15Shape.h"
#include <vector>
#include <stdexcept>

std::vector<double> RgPenta15Shape::evalH(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPenta15Shape requires NaturalCoord3d coordinates");
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

std::vector<std::vector<double>> RgPenta15Shape::evalDeriv(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPenta15Shape requires NaturalCoord3d coordinates");
    }
    
    double r = coord3d->getR();
    double s = coord3d->getS();
    double t = coord3d->getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 15 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(15));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();

    Hr[0] = -((-1 + t)*(-2 + 4 * r + 4 * s + t)) / 2.;
    Hr[1] = (-2 - 4 * r*(-1 + t) + t + t*t) / 2.;
    Hr[2] = 0;
    Hr[3] = ((-2 + 4 * r + 4 * s - t)*(1 + t)) / 2.;
    Hr[4] = ((1 + t)*(-2 + 4 * r + t)) / 2.;
    Hr[5] = 0;
    Hr[6] = 2 * (-1 + 2 * r + s)*(-1 + t);
    Hr[7] = -2 * s*(-1 + t);
    Hr[8] = 2 * s*(-1 + t);
    Hr[9] = -2 * (-1 + 2 * r + s)*(1 + t);
    Hr[10] = 2 * s*(1 + t);
    Hr[11] = -2 * s*(1 + t);
    Hr[12] = -1 + t*t;
    Hr[13] = 1 - t*t;
    Hr[14] = 0;

    Hs[0] = -((-1 + t)*(-2 + 4 * r + 4 * s + t)) / 2.;
    Hs[1] = 0;
    Hs[2] = (-2 - 4 * s*(-1 + t) + t + t*t) / 2.;
    Hs[3] = ((-2 + 4 * r + 4 * s - t)*(1 + t)) / 2.;
    Hs[4] = 0;
    Hs[5] = ((1 + t)*(-2 + 4 * s + t)) / 2.;
    Hs[6] = 2 * r*(-1 + t);
    Hs[7] = -2 * r*(-1 + t);
    Hs[8] = 2 * (-1 + r + 2 * s)*(-1 + t);
    Hs[9] = -2 * r*(1 + t);
    Hs[10] = 2 * r*(1 + t);
    Hs[11] = -2 * (-1 + r + 2 * s)*(1 + t);
    Hs[12] = -1 + t*t;
    Hs[13] = 0;
    Hs[14] = 1 - t*t;

    Ht[0] = -((-1 + r + s)*(-1 + 2 * r + 2 * s + 2 * t)) / 2.;
    Ht[1] = (r*(1 - 2 * r + 2 * t)) / 2.;
    Ht[2] = (s*(1 - 2 * s + 2 * t)) / 2.;
    Ht[3] = ((-1 + r + s)*(-1 + 2 * r + 2 * s - 2 * t)) / 2.;
    Ht[4] = r*(-0.5 + r + t);
    Ht[5] = s*(-0.5 + s + t);
    Ht[6] = 2 * r*(-1 + r + s);
    Ht[7] = -2 * r*s;
    Ht[8] = 2 * s*(-1 + r + s);
    Ht[9] = -2 * r*(-1 + r + s);
    Ht[10] = 2 * r*s;
    Ht[11] = -2 * s*(-1 + r + s);
    Ht[12] = 2 * (-1 + r + s)*t;
    Ht[13] = -2 * r*t;
    Ht[14] = -2 * s*t;
    
    return deriv;
}

std::vector<std::vector<double>> RgPenta15Shape::evalDeriv2(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPenta15Shape requires NaturalCoord3d coordinates");
    }
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 15 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(15));
    
    double* Hrr = deriv2[0].data();
    double* Hss = deriv2[1].data();
    double* Htt = deriv2[2].data();
    double* Hrs = deriv2[3].data();
    double* Hst = deriv2[4].data();
    double* Hrt = deriv2[5].data();

    Hrr[0] = 2 - 2 * t;
    Hrr[1] = 2 - 2 * t;
    Hrr[2] = 0;
    Hrr[3] = 2 * (1 + t);
    Hrr[4] = 2 * (1 + t);
    Hrr[5] = 0;
    Hrr[6] = 4 * (-1 + t);
    Hrr[7] = 0;
    Hrr[8] = 0;
    Hrr[9] = -4 * (1 + t);
    Hrr[10] = 0;
    Hrr[11] = 0;
    Hrr[12] = 0;
    Hrr[13] = 0;
    Hrr[14] = 0;

    Hss[0] = 2 - 2 * t;
    Hss[1] = 0;
    Hss[2] = 2 - 2 * t;
    Hss[3] = 2 * (1 + t);
    Hss[4] = 0;
    Hss[5] = 2 * (1 + t);
    Hss[6] = 0;
    Hss[7] = 0;
    Hss[8] = 4 * (-1 + t);
    Hss[9] = 0;
    Hss[10] = 0;
    Hss[11] = -4 * (1 + t);
    Hss[12] = 0;
    Hss[13] = 0;
    Hss[14] = 0;

    Htt[0] = 1 - r - s;
    Htt[1] = r;
    Htt[2] = s;
    Htt[3] = 1 - r - s;
    Htt[4] = r;
    Htt[5] = s;
    Htt[6] = 0;
    Htt[7] = 0;
    Htt[8] = 0;
    Htt[9] = 0;
    Htt[10] = 0;
    Htt[11] = 0;
    Htt[12] = 2 * (-1 + r + s);
    Htt[13] = -2 * r;
    Htt[14] = -2 * s;

    Hrs[0] = 2 - 2 * t;
    Hrs[1] = 0;
    Hrs[2] = 0;
    Hrs[3] = 2 * (1 + t);
    Hrs[4] = 0;
    Hrs[5] = 0;
    Hrs[6] = 2 * (-1 + t);
    Hrs[7] = 2 - 2 * t;
    Hrs[8] = 2 * (-1 + t);
    Hrs[9] = -2 * (1 + t);
    Hrs[10] = 2 * (1 + t);
    Hrs[11] = -2 * (1 + t);
    Hrs[12] = 0;
    Hrs[13] = 0;
    Hrs[14] = 0;

    Hst[0] = 1.5 - 2 * r - 2 * s - t;
    Hst[1] = 0;
    Hst[2] = 0.5 - 2 * s + t;
    Hst[3] = -1.5 + 2 * r + 2 * s - t;
    Hst[4] = 0;
    Hst[5] = -0.5 + 2 * s + t;
    Hst[6] = 2 * r;
    Hst[7] = -2 * r;
    Hst[8] = 2 * (-1 + r + 2 * s);
    Hst[9] = -2 * r;
    Hst[10] = 2 * r;
    Hst[11] = -2 * (-1 + r + 2 * s);
    Hst[12] = 2 * t;
    Hst[13] = 0;
    Hst[14] = -2 * t;

    Hrt[0] = 1.5 - 2 * r - 2 * s - t;
    Hrt[1] = 0.5 - 2 * r + t;
    Hrt[2] = 0;
    Hrt[3] = -1.5 + 2 * r + 2 * s - t;
    Hrt[4] = -0.5 + 2 * r + t;
    Hrt[5] = 0;
    Hrt[6] = 2 * (-1 + 2 * r + s);
    Hrt[7] = -2 * s;
    Hrt[8] = 2 * s;
    Hrt[9] = -2 * (-1 + 2 * r + s);
    Hrt[10] = 2 * s;
    Hrt[11] = -2 * s;
    Hrt[12] = 2 * t;
    Hrt[13] = -2 * t;
    Hrt[14] = 0;
    
    return deriv2;
}