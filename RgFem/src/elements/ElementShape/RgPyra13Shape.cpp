#include "elements/ElementShape/RgPyra13Shape.h"
#include <vector>
#include <stdexcept>

std::vector<double> RgPyra13Shape::evalH(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra13Shape requires NaturalCoord3d coordinates");
    }
    
    double r = coord3d->getR();
    double s = coord3d->getS();
    double t = coord3d->getT();
    
    std::vector<double> H(13);
    
    H[5] = 0.25*(1 - r*r)*(1 - s)*(1 - t);
    H[6] = 0.25*(1 - s*s)*(1 + r)*(1 - t);
    H[7] = 0.25*(1 - r*r)*(1 + s)*(1 - t);
    H[8] = 0.25*(1 - s*s)*(1 - r)*(1 - t);
    H[9] = 0.25*(1 - t*t)*(1 - r)*(1 - s);
    H[10] = 0.25*(1 - t*t)*(1 + r)*(1 - s);
    H[11] = 0.25*(1 - t*t)*(1 + r)*(1 + s);
    H[12] = 0.25*(1 - t*t)*(1 - r)*(1 + s);
    
    H[0] = 0.125*(1 - r)*(1 - s)*(1 - t) - 0.5*(H[5] + H[8] + H[9]);
    H[1] = 0.125*(1 + r)*(1 - s)*(1 - t) - 0.5*(H[5] + H[6] + H[10]);
    H[2] = 0.125*(1 + r)*(1 + s)*(1 - t) - 0.5*(H[6] + H[7] + H[11]);
    H[3] = 0.125*(1 - r)*(1 + s)*(1 - t) - 0.5*(H[7] + H[8] + H[12]);
    H[4] = 0.5*t*(1 + t);
    
    return H;
}

std::vector<std::vector<double>> RgPyra13Shape::evalDeriv(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra13Shape requires NaturalCoord3d coordinates");
    }
    
    double r = coord3d->getR();
    double s = coord3d->getS();
    double t = coord3d->getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 13 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(13));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();

    Hr[ 0] = 0.125 + r*(0.25 + s*(-0.25 + 0.25*t) - 0.25*t) + s*s*(-0.125 + 0.125*t) +
    s*(-0.125 + 0.125*t)*t - 0.125*t*t;
    Hr[ 1] = -0.125 + r*(0.25 + s*(-0.25 + 0.25*t) - 0.25*t) + s*s*(0.125 - 0.125*t) +
    s*(0.125 - 0.125*t)*t + 0.125*t*t;
    Hr[ 2] = -0.125 + r*(0.25 + s*(0.25 - 0.25*t) - 0.25*t) + s*s*(0.125 - 0.125*t) +
    s*(-0.125 + 0.125*t)*t + 0.125*t*t;
    Hr[ 3] = 0.125 + r*(0.25 + s*(0.25 - 0.25*t) - 0.25*t) + s*s*(-0.125 + 0.125*t) +
    s*(0.125 - 0.125*t)*t - 0.125*t*t;
    Hr[ 4] = 0;
    Hr[ 5] = -0.5*r*(-1. + s)*(-1. + t);
    Hr[ 6] = 0.25*(-1. + s*s)*(-1. + t);
    Hr[ 7] = 0.5*r*(1. + s)*(-1. + t);
    Hr[ 8] = -0.25*(-1. + s*s)*(-1. + t);
    Hr[ 9] = -0.25*(-1. + s)*(-1. + t*t);
    Hr[10] = 0.25*(-1. + s)*(-1. + t*t);
    Hr[11] = -0.25*(1. + s)*(-1. + t*t);
    Hr[12] = 0.25*(1. + s)*(-1. + t*t);
    
    Hs[ 0] = 0.125 + s*(0.25 - 0.25*t) + r*r*(-0.125 + 0.125*t) - 0.125*t*t +
    r*(s*(-0.25 + 0.25*t) + (-0.125 + 0.125*t)*t);
    Hs[ 1] = 0.125 + s*(0.25 - 0.25*t) + r*r*(-0.125 + 0.125*t) - 0.125*t*t +
    r*(s*(0.25 - 0.25*t) + (0.125 - 0.125*t)*t);
    Hs[ 2] = -0.125 + s*(0.25 - 0.25*t) + r*r*(0.125 - 0.125*t) + 0.125*t*t +
    r*(s*(0.25 - 0.25*t) + (-0.125 + 0.125*t)*t);
    Hs[ 3] = -0.125 + s*(0.25 - 0.25*t) + r*r*(0.125 - 0.125*t) + 0.125*t*t +
    r*(s*(-0.25 + 0.25*t) + (0.125 - 0.125*t)*t);
    Hs[ 4] = 0;
    Hs[ 5] = -0.25*(-1. + r*r)*(-1. + t);
    Hs[ 6] = 0.5*(1. + r)*s*(-1. + t);
    Hs[ 7] = 0.25*(-1. + r*r)*(-1. + t);
    Hs[ 8] = -0.5*(-1. + r)*s*(-1. + t);
    Hs[ 9] = -0.25*(-1. + r)*(-1. + t*t);
    Hs[10] = 0.25*(1. + r)*(-1. + t*t);
    Hs[11] = -0.25*(1. + r)*(-1. + t*t);
    Hs[12] = 0.25*(-1. + r)*(-1. + t*t);

    Ht[ 0] = -0.125*(-1. + r)*(-1. + s) + 0.125*(-1. + r*r)*(-1. + s) +
    0.125*(-1. + r)*(-1. + s*s) + 0.25*(-1. + r)*(-1. + s)*t;
    Ht[ 1] = 0.125*(1. + r)*(-1. + s) + 0.125*(-1. + r*r)*(-1. + s) -
    0.125*(1. + r)*(-1. + s*s) - 0.25*(1. + r)*(-1. + s)*t;
    Ht[ 2] = -0.125*(1. + r)*(1. + s) - 0.125*(-1. + r*r)*(1. + s) -
    0.125*(1. + r)*(-1. + s*s) + 0.25*(1. + r)*(1. + s)*t;
    Ht[ 3] = 0.125*(-1. + r)*(1. + s) - 0.125*(-1. + r*r)*(1. + s) +
    0.125*(-1. + r)*(-1. + s*s) - 0.25*(-1. + r)*(1. + s)*t;
    Ht[ 4] = 0.5 + 1.*t;
    Ht[ 5] = -0.25*(-1. + r*r)*(-1. + s);
    Ht[ 6] = 0.25*(1. + r)*(-1. + s*s);
    Ht[ 7] = 0.25*(-1. + r*r)*(1. + s);
    Ht[ 8] = -0.25*(-1. + r)*(-1. + s*s);
    Ht[ 9] = -0.5*(-1. + r)*(-1. + s)*t;
    Ht[10] = 0.5*(1. + r)*(-1. + s)*t;
    Ht[11] = -0.5*(1. + r)*(1. + s)*t;
    Ht[12] = 0.5*(-1. + r)*(1. + s)*t;
    
    return deriv;
}

std::vector<std::vector<double>> RgPyra13Shape::evalDeriv2(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra13Shape requires NaturalCoord3d coordinates");
    }
    
    double r = coord3d->getR();
    double s = coord3d->getS();
    double t = coord3d->getT();
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 13 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(13));
    
    double* Hrr = deriv2[0].data();
    double* Hss = deriv2[1].data();
    double* Htt = deriv2[2].data();
    double* Hrs = deriv2[3].data();
    double* Hst = deriv2[4].data();
    double* Hrt = deriv2[5].data();

    Hrr[ 0] = 0.25*(-1. + s)*(-1. + t);
    Hrr[ 1] = 0.25*(-1. + s)*(-1. + t);
    Hrr[ 2] = -0.25*(1. + s)*(-1. + t);
    Hrr[ 3] = -0.25*(1. + s)*(-1. + t);
    Hrr[ 4] = 0;
    Hrr[ 5] = -0.5*(-1. + s)*(-1. + t);
    Hrr[ 6] = 0.;
    Hrr[ 7] = 0.5*(1. + s)*(-1. + t);
    Hrr[ 8] = 0.;
    Hrr[ 9] = 0.;
    Hrr[10] = 0.;
    Hrr[11] = 0.;
    Hrr[12] = 0.;

    Hss[ 0] = 0.25*(-1. + r)*(-1. + t);
    Hss[ 1] = -0.25*(1. + r)*(-1. + t);
    Hss[ 2] = -0.25*(1. + r)*(-1. + t);
    Hss[ 3] = 0.25*(-1. + r)*(-1. + t);
    Hss[ 4] = 0;
    Hss[ 5] = 0;
    Hss[ 6] = 0.5*(1. + r)*(-1. + t);
    Hss[ 7] = 0;
    Hss[ 8] = -0.5*(-1. + r)*(-1. + t);
    Hss[ 9] = 0;
    Hss[10] = 0;
    Hss[11] = 0;
    Hss[12] = 0;

    Htt[ 0] = 0.25*(-1. + r)*(-1. + s);
    Htt[ 1] = -0.25*(1. + r)*(-1. + s);
    Htt[ 2] = 0.25*(1. + r)*(1. + s);
    Htt[ 3] = -0.25*(-1. + r)*(1. + s);
    Htt[ 4] = 1.;
    Htt[ 5] = 0;
    Htt[ 6] = 0;
    Htt[ 7] = 0;
    Htt[ 8] = 0;
    Htt[ 9] = -0.5*(-1. + r)*(-1. + s);
    Htt[10] = 0.5*(1. + r)*(-1. + s);
    Htt[11] = -0.5*(1. + r)*(1. + s);
    Htt[12] = 0.5*(-1. + r)*(1. + s);

    Hrs[ 0] = r*(-0.25 + 0.25*t) + s*(-0.25 + 0.25*t) - 0.125*t + 0.125*t*t;
    Hrs[ 1] = s*(0.25 - 0.25*t) + r*(-0.25 + 0.25*t) + 0.125*t - 0.125*t*t;
    Hrs[ 2] = r*(0.25 - 0.25*t) + s*(0.25 - 0.25*t) - 0.125*t + 0.125*t*t;
    Hrs[ 3] = r*(0.25 - 0.25*t) + s*(-0.25 + 0.25*t) + 0.125*t - 0.125*t*t;
    Hrs[ 4] = 0;
    Hrs[ 5] = -0.5*r*(-1. + t);
    Hrs[ 6] = 0.5*s*(-1. + t);
    Hrs[ 7] = 0.5*r*(-1. + t);
    Hrs[ 8] = -0.5*s*(-1. + t);
    Hrs[ 9] = -0.25*(-1. + t*t);
    Hrs[10] = 0.25*(-1. + t*t);
    Hrs[11] = -0.25*(-1. + t*t);
    Hrs[12] = 0.25*(-1. + t*t);
    
    Hrs[ 0] = 0.125*r*r - 0.25*s + r*(-0.125 + 0.25*s + 0.25*t) - 0.25*t;
    Hst[ 1] = 0.125*r*r - 0.25*s + r*(0.125 - 0.25*s - 0.25*t) - 0.25*t;
    Hst[ 2] = -0.125*r*r - 0.25*s + r*(-0.125 - 0.25*s + 0.25*t) + 0.25*t;
    Hst[ 3] = -0.125*r*r - 0.25*s + r*(0.125 + 0.25*s - 0.25*t) + 0.25*t;
    Hst[ 4] = 0;
    Hst[ 5] = -0.25*(-1. + r*r);
    Hst[ 6] = 0.5*(1. + r)*s;
    Hst[ 7] = 0.25*(-1. + r*r);
    Hst[ 8] = -0.5*(-1. + r)*s;
    Hst[ 9] = -0.5*(-1. + r)*t;
    Hst[10] = 0.5*(1. + r)*t;
    Hst[11] = -0.5*(1. + r)*t;
    Hst[12] = 0.5*(-1. + r)*t;
    
    Hrt[ 0] = r*(-0.25 + 0.25*s) + 0.125*s*s + s*(-0.125 + 0.25*t) - 0.25*t;
    Hrt[ 1] = r*(-0.25 + 0.25*s) - 0.125*s*s + s*(0.125 - 0.25*t) + 0.25*t;
    Hrt[ 2] = r*(-0.25 - 0.25*s) - 0.125*s*s + s*(-0.125 + 0.25*t) + 0.25*t;
    Hrt[ 3] = r*(-0.25 - 0.25*s) + 0.125*s*s + s*(0.125 - 0.25*t) - 0.25*t;
    Hrt[ 4] = 0;
    Hrt[ 5] = -0.5*r*(-1. + s);
    Hrt[ 6] = 0.25*(-1. + s*s);
    Hrt[ 7] = 0.5*r*(1. + s);
    Hrt[ 8] = -0.25*(-1. + s*s);
    Hrt[ 9] = -0.5*(-1. + s)*t;
    Hrt[10] = 0.5*(-1. + s)*t;
    Hrt[11] = -0.5*(1. + s)*t;
    Hrt[12] = 0.5*(1. + s)*t;

    Hrs[0] = 0.125*(1.0 - t); Hrt[0] = 0.125*(1.0 - s); Hst[0] = 0.125*(1.0 - r);
    Hrs[1] = -0.125*(1.0 - t); Hrt[1] = -0.125*(1.0 - s); Hst[1] = 0.125*(1.0 + r);
    Hrs[2] = 0.125*(1.0 - t); Hrt[2] = -0.125*(1.0 + s); Hst[2] = -0.125*(1.0 + r);
    Hrs[3] = -0.125*(1.0 - t); Hrt[3] = 0.125*(1.0 + s); Hst[3] = -0.125*(1.0 - r);
    Hrs[4] = 0.0; Hrt[4] = 0.0; Hst[4] = 0.0;
    
    return deriv2;
}