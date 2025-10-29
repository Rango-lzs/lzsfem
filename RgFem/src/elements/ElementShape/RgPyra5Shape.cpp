#include "elements/ElementShape/RgPyra5Shape.h"
#include <vector>
#include <stdexcept>

std::vector<double> RgPyra5Shape::evalH(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra5Shape requires NaturalCoord3d coordinates");
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

std::vector<std::vector<double>> RgPyra5Shape::evalDeriv(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra5Shape requires NaturalCoord3d coordinates");
    }
    
    double r = coord3d->getR();
    double s = coord3d->getS();
    double t = coord3d->getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 5 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(5));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();

    Hr[0] = -0.125*(1.0 - s)*(1.0 - t);
    Hr[1] = 0.125*(1.0 - s)*(1.0 - t);
    Hr[2] = 0.125*(1.0 + s)*(1.0 - t);
    Hr[3] = -0.125*(1.0 + s)*(1.0 - t);
    Hr[4] = 0.0;

    Hs[0] = -0.125*(1.0 - r)*(1.0 - t);
    Hs[1] = -0.125*(1.0 + r)*(1.0 - t);
    Hs[2] = 0.125*(1.0 + r)*(1.0 - t);
    Hs[3] = 0.125*(1.0 - r)*(1.0 - t);
    Hs[4] = 0.0;

    Ht[0] = -0.125*(1.0 - r)*(1.0 - s);
    Ht[1] = -0.125*(1.0 + r)*(1.0 - s);
    Ht[2] = -0.125*(1.0 + r)*(1.0 + s);
    Ht[3] = -0.125*(1.0 - r)*(1.0 + s);
    Ht[4] = 0.5;
    
    return deriv;
}

std::vector<std::vector<double>> RgPyra5Shape::evalDeriv2(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgPyra5Shape requires NaturalCoord3d coordinates");
    }
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 5 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(5, 0.0));
    
    double* Hrr = deriv2[0].data();
    double* Hss = deriv2[1].data();
    double* Htt = deriv2[2].data();
    double* Hrs = deriv2[3].data();
    double* Hst = deriv2[4].data();
    double* Hrt = deriv2[5].data();

    Hrr[0] = 0.0; Hss[0] = 0.0; Htt[0] = 0.0;
    Hrr[1] = 0.0; Hss[1] = 0.0; Htt[1] = 0.0;
    Hrr[2] = 0.0; Hss[2] = 0.0; Htt[2] = 0.0;
    Hrr[3] = 0.0; Hss[3] = 0.0; Htt[3] = 0.0;
    Hrr[4] = 0.0; Hss[4] = 0.0; Htt[4] = 0.0;

    Hrs[0] = 0.125*(1.0 - t); Hrt[0] = 0.125*(1.0 - s); Hst[0] = 0.125*(1.0 - r);
    Hrs[1] = -0.125*(1.0 - t); Hrt[1] = -0.125*(1.0 - s); Hst[1] = 0.125*(1.0 + r);
    Hrs[2] = 0.125*(1.0 - t); Hrt[2] = -0.125*(1.0 + s); Hst[2] = -0.125*(1.0 + r);
    Hrs[3] = -0.125*(1.0 - t); Hrt[3] = 0.125*(1.0 + s); Hst[3] = -0.125*(1.0 - r);
    Hrs[4] = 0.0; Hrt[4] = 0.0; Hst[4] = 0.0;
    
    return deriv2;
}