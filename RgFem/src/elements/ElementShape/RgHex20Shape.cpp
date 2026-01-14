#include "elements/ElementShape/RgHex20Shape.h"
#include <vector>
#include <stdexcept>
#include "../NaturalCoord.h"


std::vector<double> RgHex20Shape::evalH(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    
    std::vector<double> H(20);
    
    // Midside nodes
    H[8] = 0.25*(1 - r*r)*(1 - s)*(1 - t);
    H[9] = 0.25*(1 - s*s)*(1 + r)*(1 - t);
    H[10] = 0.25*(1 - r*r)*(1 + s)*(1 - t);
    H[11] = 0.25*(1 - s*s)*(1 - r)*(1 - t);
    H[12] = 0.25*(1 - r*r)*(1 - s)*(1 + t);
    H[13] = 0.25*(1 - s*s)*(1 + r)*(1 + t);
    H[14] = 0.25*(1 - r*r)*(1 + s)*(1 + t);
    H[15] = 0.25*(1 - s*s)*(1 - r)*(1 + t);
    H[16] = 0.25*(1 - t*t)*(1 - r)*(1 - s);
    H[17] = 0.25*(1 - t*t)*(1 + r)*(1 - s);
    H[18] = 0.25*(1 - t*t)*(1 + r)*(1 + s);
    H[19] = 0.25*(1 - t*t)*(1 - r)*(1 + s);

    // Corner nodes
    H[0] = 0.125*(1 - r)*(1 - s)*(1 - t) - 0.5*(H[8] + H[11] + H[16]);
    H[1] = 0.125*(1 + r)*(1 - s)*(1 - t) - 0.5*(H[8] + H[9] + H[17]);
    H[2] = 0.125*(1 + r)*(1 + s)*(1 - t) - 0.5*(H[9] + H[10] + H[18]);
    H[3] = 0.125*(1 - r)*(1 + s)*(1 - t) - 0.5*(H[10] + H[11] + H[19]);
    H[4] = 0.125*(1 - r)*(1 - s)*(1 + t) - 0.5*(H[12] + H[15] + H[16]);
    H[5] = 0.125*(1 + r)*(1 - s)*(1 + t) - 0.5*(H[12] + H[13] + H[17]);
    H[6] = 0.125*(1 + r)*(1 + s)*(1 + t) - 0.5*(H[13] + H[14] + H[18]);
    H[7] = 0.125*(1 - r)*(1 + s)*(1 + t) - 0.5*(H[14] + H[15] + H[19]);
    
    return H;
}

std::vector<std::vector<double>> RgHex20Shape::evalDeriv(const NaturalCoord& coord)
{   
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 20 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(20));
    
    double* Hr = deriv[0].data();
    double* Hs = deriv[1].data();
    double* Ht = deriv[2].data();
    
    // dH/dr
    Hr[8] = -0.5*r*(1 - s)*(1 - t);
    Hr[9] = 0.25*(1 - s*s)*(1 - t);
    Hr[10] = -0.5*r*(1 + s)*(1 - t);
    Hr[11] = -0.25*(1 - s*s)*(1 - t);
    Hr[12] = -0.5*r*(1 - s)*(1 + t);
    Hr[13] = 0.25*(1 - s*s)*(1 + t);
    Hr[14] = -0.5*r*(1 + s)*(1 + t);
    Hr[15] = -0.25*(1 - s*s)*(1 + t);
    Hr[16] = -0.25*(1 - t*t)*(1 - s);
    Hr[17] = 0.25*(1 - t*t)*(1 - s);
    Hr[18] = 0.25*(1 - t*t)*(1 + s);
    Hr[19] = -0.25*(1 - t*t)*(1 + s);

    Hr[0] = -0.125*(1 - s)*(1 - t) - 0.5*(Hr[8] + Hr[11] + Hr[16]);
    Hr[1] = 0.125*(1 - s)*(1 - t) - 0.5*(Hr[8] + Hr[9] + Hr[17]);
    Hr[2] = 0.125*(1 + s)*(1 - t) - 0.5*(Hr[9] + Hr[10] + Hr[18]);
    Hr[3] = -0.125*(1 + s)*(1 - t) - 0.5*(Hr[10] + Hr[11] + Hr[19]);
    Hr[4] = -0.125*(1 - s)*(1 + t) - 0.5*(Hr[12] + Hr[15] + Hr[16]);
    Hr[5] = 0.125*(1 - s)*(1 + t) - 0.5*(Hr[12] + Hr[13] + Hr[17]);
    Hr[6] = 0.125*(1 + s)*(1 + t) - 0.5*(Hr[13] + Hr[14] + Hr[18]);
    Hr[7] = -0.125*(1 + s)*(1 + t) - 0.5*(Hr[14] + Hr[15] + Hr[19]);

    // dH/ds
    Hs[8] = -0.25*(1 - r*r)*(1 - t);
    Hs[9] = -0.5*s*(1 + r)*(1 - t);
    Hs[10] = 0.25*(1 - r*r)*(1 - t);
    Hs[11] = -0.5*s*(1 - r)*(1 - t);
    Hs[12] = -0.25*(1 - r*r)*(1 + t);
    Hs[13] = -0.5*s*(1 + r)*(1 + t);
    Hs[14] = 0.25*(1 - r*r)*(1 + t);
    Hs[15] = -0.5*s*(1 - r)*(1 + t);
    Hs[16] = -0.25*(1 - t*t)*(1 - r);
    Hs[17] = -0.25*(1 - t*t)*(1 + r);
    Hs[18] = 0.25*(1 - t*t)*(1 + r);
    Hs[19] = 0.25*(1 - t*t)*(1 - r);

    Hs[0] = -0.125*(1 - r)*(1 - t) - 0.5*(Hs[8] + Hs[11] + Hs[16]);
    Hs[1] = -0.125*(1 + r)*(1 - t) - 0.5*(Hs[8] + Hs[9] + Hs[17]);
    Hs[2] = 0.125*(1 + r)*(1 - t) - 0.5*(Hs[9] + Hs[10] + Hs[18]);
    Hs[3] = 0.125*(1 - r)*(1 - t) - 0.5*(Hs[10] + Hs[11] + Hs[19]);
    Hs[4] = -0.125*(1 - r)*(1 + t) - 0.5*(Hs[12] + Hs[15] + Hs[16]);
    Hs[5] = -0.125*(1 + r)*(1 + t) - 0.5*(Hs[12] + Hs[13] + Hs[17]);
    Hs[6] = 0.125*(1 + r)*(1 + t) - 0.5*(Hs[13] + Hs[14] + Hs[18]);
    Hs[7] = 0.125*(1 - r)*(1 + t) - 0.5*(Hs[14] + Hs[15] + Hs[19]);

    // dH/dt
    Ht[8] = -0.25*(1 - r*r)*(1 - s);
    Ht[9] = -0.25*(1 - s*s)*(1 + r);
    Ht[10] = -0.25*(1 - r*r)*(1 + s);
    Ht[11] = -0.25*(1 - s*s)*(1 - r);
    Ht[12] = 0.25*(1 - r*r)*(1 - s);
    Ht[13] = 0.25*(1 - s*s)*(1 + r);
    Ht[14] = 0.25*(1 - r*r)*(1 + s);
    Ht[15] = 0.25*(1 - s*s)*(1 - r);
    Ht[16] = -0.5*t*(1 - r)*(1 - s);
    Ht[17] = -0.5*t*(1 + r)*(1 - s);
    Ht[18] = -0.5*t*(1 + r)*(1 + s);
    Ht[19] = -0.5*t*(1 - r)*(1 + s);

    Ht[0] = -0.125*(1 - r)*(1 - s) - 0.5*(Ht[8] + Ht[11] + Ht[16]);
    Ht[1] = -0.125*(1 + r)*(1 - s) - 0.5*(Ht[8] + Ht[9] + Ht[17]);
    Ht[2] = -0.125*(1 + r)*(1 + s) - 0.5*(Ht[9] + Ht[10] + Ht[18]);
    Ht[3] = -0.125*(1 - r)*(1 + s) - 0.5*(Ht[10] + Ht[11] + Ht[19]);
    Ht[4] = 0.125*(1 - r)*(1 - s) - 0.5*(Ht[12] + Ht[15] + Ht[16]);
    Ht[5] = 0.125*(1 + r)*(1 - s) - 0.5*(Ht[12] + Ht[13] + Ht[17]);
    Ht[6] = 0.125*(1 + r)*(1 + s) - 0.5*(Ht[13] + Ht[14] + Ht[18]);
    Ht[7] = 0.125*(1 - r)*(1 + s) - 0.5*(Ht[14] + Ht[15] + Ht[19]);
    
    return deriv;
}

std::vector<std::vector<double>> RgHex20Shape::evalDeriv2(const NaturalCoord& coord)
{
    double r = coord.getR();
    double s = coord.getS();
    double t = coord.getT();
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 20 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(20));
    
    double* Hrr = deriv2[0].data();
    double* Hss = deriv2[1].data();
    double* Htt = deriv2[2].data();
    double* Hrs = deriv2[3].data();
    double* Hst = deriv2[4].data();
    double* Hrt = deriv2[5].data();
    
    // Hrr
    Hrr[8] = -0.5*(1 - s)*(1 - t);
    Hrr[9] = 0.;
    Hrr[10] = -0.5*(1 + s)*(1 - t);
    Hrr[11] = 0.;
    Hrr[12] = -0.5*(1 - s)*(1 + t);
    Hrr[13] = 0.;
    Hrr[14] = -0.5*(1 + s)*(1 + t);
    Hrr[15] = 0.;
    Hrr[16] = 0.;
    Hrr[17] = 0.;
    Hrr[18] = 0.;
    Hrr[19] = 0.;

    Hrr[0] = -0.5*(Hrr[8] + Hrr[11] + Hrr[16]);
    Hrr[1] = -0.5*(Hrr[8] + Hrr[9] + Hrr[17]);
    Hrr[2] = -0.5*(Hrr[9] + Hrr[10] + Hrr[18]);
    Hrr[3] = -0.5*(Hrr[10] + Hrr[11] + Hrr[19]);
    Hrr[4] = -0.5*(Hrr[12] + Hrr[15] + Hrr[16]);
    Hrr[5] = -0.5*(Hrr[12] + Hrr[13] + Hrr[17]);
    Hrr[6] = -0.5*(Hrr[13] + Hrr[14] + Hrr[18]);
    Hrr[7] = -0.5*(Hrr[14] + Hrr[15] + Hrr[19]);

    // Hss
    Hss[8] = 0.;
    Hss[9] = -0.5*(1 + r)*(1 - t);
    Hss[10] = 0.;
    Hss[11] = -0.5*(1 - r)*(1 - t);
    Hss[12] = 0.;
    Hss[13] = -0.5*(1 + r)*(1 + t);
    Hss[14] = 0.;
    Hss[15] = -0.5*(1 - r)*(1 + t);
    Hss[16] = 0.;
    Hss[17] = 0.;
    Hss[18] = 0.;
    Hss[19] = 0.;

    Hss[0] = -0.5*(Hss[8] + Hss[11] + Hss[16]);
    Hss[1] = -0.5*(Hss[8] + Hss[9] + Hss[17]);
    Hss[2] = -0.5*(Hss[9] + Hss[10] + Hss[18]);
    Hss[3] = -0.5*(Hss[10] + Hss[11] + Hss[19]);
    Hss[4] = -0.5*(Hss[12] + Hss[15] + Hss[16]);
    Hss[5] = -0.5*(Hss[12] + Hss[13] + Hss[17]);
    Hss[6] = -0.5*(Hss[13] + Hss[14] + Hss[18]);
    Hss[7] = -0.5*(Hss[14] + Hss[15] + Hss[19]);

    // Htt
    Htt[8] = 0.;
    Htt[9] = 0.;
    Htt[10] = 0.;
    Htt[11] = 0.;
    Htt[12] = 0.;
    Htt[13] = 0.;
    Htt[14] = 0.;
    Htt[15] = 0.;
    Htt[16] = -0.5*(1 - r)*(1 - s);
    Htt[17] = -0.5*(1 + r)*(1 - s);
    Htt[18] = -0.5*(1 + r)*(1 + s);
    Htt[19] = -0.5*(1 - r)*(1 + s);

    Htt[0] = -0.5*(Htt[8] + Htt[11] + Htt[16]);
    Htt[1] = -0.5*(Htt[8] + Htt[9] + Htt[17]);
    Htt[2] = -0.5*(Htt[9] + Htt[10] + Htt[18]);
    Htt[3] = -0.5*(Htt[10] + Htt[11] + Htt[19]);
    Htt[4] = -0.5*(Htt[12] + Htt[15] + Htt[16]);
    Htt[5] = -0.5*(Htt[12] + Htt[13] + Htt[17]);
    Htt[6] = -0.5*(Htt[13] + Htt[14] + Htt[18]);
    Htt[7] = -0.5*(Htt[14] + Htt[15] + Htt[19]);

    // Hrs
    Hrs[8] = 0.5*r*(1 - t);
    Hrs[9] = -0.5*s*(1 - t);
    Hrs[10] = -0.5*r*(1 - t);
    Hrs[11] = 0.5*s*(1 - t);
    Hrs[12] = 0.5*r*(1 + t);
    Hrs[13] = -0.5*s*(1 + t);
    Hrs[14] = -0.5*r*(1 + t);
    Hrs[15] = 0.5*s*(1 + t);
    Hrs[16] = 0.25*(1 - t*t);
    Hrs[17] = -0.25*(1 - t*t);
    Hrs[18] = 0.25*(1 - t*t);
    Hrs[19] = -0.25*(1 - t*t);

    Hrs[0] = 0.125*(1 - t) - 0.5*(Hrs[8] + Hrs[11] + Hrs[16]);
    Hrs[1] = -0.125*(1 - t) - 0.5*(Hrs[8] + Hrs[9] + Hrs[17]);
    Hrs[2] = 0.125*(1 - t) - 0.5*(Hrs[9] + Hrs[10] + Hrs[18]);
    Hrs[3] = -0.125*(1 - t) - 0.5*(Hrs[10] + Hrs[11] + Hrs[19]);
    Hrs[4] = 0.125*(1 + t) - 0.5*(Hrs[12] + Hrs[15] + Hrs[16]);
    Hrs[5] = -0.125*(1 + t) - 0.5*(Hrs[12] + Hrs[13] + Hrs[17]);
    Hrs[6] = 0.125*(1 + t) - 0.5*(Hrs[13] + Hrs[14] + Hrs[18]);
    Hrs[7] = -0.125*(1 + t) - 0.5*(Hrs[14] + Hrs[15] + Hrs[19]);

    // Hst
    Hst[8] = 0.25*(1 - r*r);
    Hst[9] = 0.5*s*(1 + r);
    Hst[10] = -0.25*(1 - r*r);
    Hst[11] = 0.5*s*(1 - r);
    Hst[12] = -0.25*(1 - r*r);
    Hst[13] = -0.5*s*(1 + r);
    Hst[14] = 0.25*(1 - r*r);
    Hst[15] = -0.5*s*(1 - r);
    Hst[16] = 0.5*t*(1 - r);
    Hst[17] = 0.5*t*(1 + r);
    Hst[18] = -0.5*t*(1 + r);
    Hst[19] = -0.5*t*(1 - r);

    Hst[0] = 0.125*(1 - r) - 0.5*(Hst[8] + Hst[11] + Hst[16]);
    Hst[1] = 0.125*(1 + r) - 0.5*(Hst[8] + Hst[9] + Hst[17]);
    Hst[2] = -0.125*(1 + r) - 0.5*(Hst[9] + Hst[10] + Hst[18]);
    Hst[3] = -0.125*(1 - r) - 0.5*(Hst[10] + Hst[11] + Hst[19]);
    Hst[4] = -0.125*(1 - r) - 0.5*(Hst[12] + Hst[15] + Hst[16]);
    Hst[5] = -0.125*(1 + r) - 0.5*(Hst[12] + Hst[13] + Hst[17]);
    Hst[6] = 0.125*(1 + r) - 0.5*(Hst[13] + Hst[14] + Hst[18]);
    Hst[7] = 0.125*(1 - r) - 0.5*(Hst[14] + Hst[15] + Hst[19]);

    // Hrt
    Hrt[8] = 0.5*r*(1 - s);
    Hrt[9] = -0.25*(1 - s*s);
    Hrt[10] = 0.5*r*(1 + s);
    Hrt[11] = 0.25*(1 - s*s);
    Hrt[12] = -0.5*r*(1 - s);
    Hrt[13] = 0.25*(1 - s*s);
    Hrt[14] = -0.5*r*(1 + s);
    Hrt[15] = -0.25*(1 - s*s);
    Hrt[16] = 0.5*t*(1 - s);
    Hrt[17] = -0.5*t*(1 - s);
    Hrt[18] = -0.5*t*(1 + s);
    Hrt[19] = 0.5*t*(1 + s);

    Hrt[0] = 0.125*(1 - s) - 0.5*(Hrt[8] + Hrt[11] + Hrt[16]);
    Hrt[1] = -0.125*(1 - s) - 0.5*(Hrt[8] + Hrt[9] + Hrt[17]);
    Hrt[2] = -0.125*(1 + s) - 0.5*(Hrt[9] + Hrt[10] + Hrt[18]);
    Hrt[3] = 0.125*(1 + s) - 0.5*(Hrt[10] + Hrt[11] + Hrt[19]);
    Hrt[4] = -0.125*(1 - s) - 0.5*(Hrt[12] + Hrt[15] + Hrt[16]);
    Hrt[5] = 0.125*(1 - s) - 0.5*(Hrt[12] + Hrt[13] + Hrt[17]);
    Hrt[6] = 0.125*(1 + s) - 0.5*(Hrt[13] + Hrt[14] + Hrt[18]);
    Hrt[7] = -0.125*(1 + s) - 0.5*(Hrt[14] + Hrt[15] + Hrt[19]);
    
    return deriv2;
}