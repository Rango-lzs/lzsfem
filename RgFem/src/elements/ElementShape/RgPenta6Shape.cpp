#include "elements/ElementShape/RgPenta6Shape.h"
#include <vector>
#include "../NaturalCoord.h"



std::vector<double> RgPenta6Shape::evalH(const NaturalCoord& coord)
{
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    
    std::vector<double> H(6);
    
    H[0] = 0.5*(1 - t)*(1 - r - s);
    H[1] = 0.5*(1 - t)*r;
    H[2] = 0.5*(1 - t)*s;
    H[3] = 0.5*(1 + t)*(1 - r - s);
    H[4] = 0.5*(1 + t)*r;
    H[5] = 0.5*(1 + t)*s;
    
    return H;
}

std::vector<std::vector<double>> RgPenta6Shape::evalDeriv(const NaturalCoord& coord)
{
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 6 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(6));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();

    Hr[0] = -0.5*(1 - t);
    Hr[1] = 0.5*(1 - t);
    Hr[2] = 0.0;
    Hr[3] = -0.5*(1 + t);
    Hr[4] = 0.5*(1 + t);
    Hr[5] = 0.0;

    Hs[0] = -0.5*(1 - t);
    Hs[1] = 0.0;
    Hs[2] = 0.5*(1 - t);
    Hs[3] = -0.5*(1 + t);
    Hs[4] = 0.0;
    Hs[5] = 0.5*(1 + t);

    Ht[0] = -0.5*(1 - r - s);
    Ht[1] = -0.5*r;
    Ht[2] = -0.5*s;
    Ht[3] = 0.5*(1 - r - s);
    Ht[4] = 0.5*r;
    Ht[5] = 0.5*s;
    
    return deriv;
}

std::vector<std::vector<double>> RgPenta6Shape::evalDeriv2(const NaturalCoord& coord)
{
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 6 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(6, 0.0));
    
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
    Hrr[5] = 0.0; Hss[5] = 0.0; Htt[5] = 0.0;

    Hrs[0] = 0.0; Hst[0] = 0.5; Hrt[0] = 0.5;
    Hrs[1] = 0.0; Hst[1] = 0.0; Hrt[1] = -0.5;
    Hrs[2] = 0.0; Hst[2] = -0.5; Hrt[2] = 0.0;
    Hrs[3] = 0.0; Hst[3] = -0.5; Hrt[3] = -0.5;
    Hrs[4] = 0.0; Hst[4] = 0.0; Hrt[4] = 0.5;
    Hrs[5] = 0.0; Hst[5] = 0.5; Hrt[5] = 0.0;
    
    return deriv2;
}