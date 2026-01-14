#include "elements/ElementShape/RgTet10Shape.h"
#include <vector>
#include "../NaturalCoord.h"



std::vector<double> RgTet10Shape::evalH(const NaturalCoord& coord)
{
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    double u = 1.0 - r - s - t;
    
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

std::vector<std::vector<double>> RgTet10Shape::evalDeriv(const NaturalCoord& coord)
{
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    double u = 1.0 - r - s - t;
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 10 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(10));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();

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
    
    return deriv;
}

std::vector<std::vector<double>> RgTet10Shape::evalDeriv2(const NaturalCoord& coord)
{
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 10 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(10));
    
    double* Hrr = deriv2[0].data();
    double* Hss = deriv2[1].data();
    double* Htt = deriv2[2].data();
    double* Hrs = deriv2[3].data();
    double* Hst = deriv2[4].data();
    double* Hrt = deriv2[5].data();

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
    
    return deriv2;
}