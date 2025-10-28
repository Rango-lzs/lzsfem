#include "elements/ElementShape/RgHex8Shape.h"
#include <cassert>

std::vector<double> RgHex8shape::evalH(const NaturalCoord& coord)
{
    // Cast to NaturalCoord3d to access r, s, t coordinates
    const NaturalCoord3d& coord3d = static_cast<const NaturalCoord3d&>(coord);
    double r = coord3d.getR();
    double s = coord3d.getS();
    double t = coord3d.getT();
    
    std::vector<double> H(8);
    
    H[0] = 0.125*(1 - r)*(1 - s)*(1 - t);
    H[1] = 0.125*(1 + r)*(1 - s)*(1 - t);
    H[2] = 0.125*(1 + r)*(1 + s)*(1 - t);
    H[3] = 0.125*(1 - r)*(1 + s)*(1 - t);
    H[4] = 0.125*(1 - r)*(1 - s)*(1 + t);
    H[5] = 0.125*(1 + r)*(1 - s)*(1 + t);
    H[6] = 0.125*(1 + r)*(1 + s)*(1 + t);
    H[7] = 0.125*(1 - r)*(1 + s)*(1 + t);
    
    return H;
}

std::vector<std::vector<double>> RgHex8shape::evalDeriv(const NaturalCoord& coord)
{
    const NaturalCoord3d& coord3d = static_cast<const NaturalCoord3d&>(coord);
    double r = coord3d.getR();
    double s = coord3d.getS();
    double t = coord3d.getT();
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 8 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(8));
    
    // dH/dr
    deriv[0][0] = -0.125*(1 - s)*(1 - t);
    deriv[0][1] =  0.125*(1 - s)*(1 - t);
    deriv[0][2] =  0.125*(1 + s)*(1 - t);
    deriv[0][3] = -0.125*(1 + s)*(1 - t);
    deriv[0][4] = -0.125*(1 - s)*(1 + t);
    deriv[0][5] =  0.125*(1 - s)*(1 + t);
    deriv[0][6] =  0.125*(1 + s)*(1 + t);
    deriv[0][7] = -0.125*(1 + s)*(1 + t);
    
    // dH/ds
    deriv[1][0] = -0.125*(1 - r)*(1 - t);
    deriv[1][1] = -0.125*(1 + r)*(1 - t);
    deriv[1][2] =  0.125*(1 + r)*(1 - t);
    deriv[1][3] =  0.125*(1 - r)*(1 - t);
    deriv[1][4] = -0.125*(1 - r)*(1 + t);
    deriv[1][5] = -0.125*(1 + r)*(1 + t);
    deriv[1][6] =  0.125*(1 + r)*(1 + t);
    deriv[1][7] =  0.125*(1 - r)*(1 + t);
    
    // dH/dt
    deriv[2][0] = -0.125*(1 - r)*(1 - s);
    deriv[2][1] = -0.125*(1 + r)*(1 - s);
    deriv[2][2] = -0.125*(1 + r)*(1 + s);
    deriv[2][3] = -0.125*(1 - r)*(1 + s);
    deriv[2][4] =  0.125*(1 - r)*(1 - s);
    deriv[2][5] =  0.125*(1 + r)*(1 - s);
    deriv[2][6] =  0.125*(1 + r)*(1 + s);
    deriv[2][7] =  0.125*(1 - r)*(1 + s);
    
    return deriv;
}

std::vector<std::vector<double>> RgHex8shape::evalDeriv2(const NaturalCoord& coord)
{
    const NaturalCoord3d& coord3d = static_cast<const NaturalCoord3d&>(coord);
    double r = coord3d.getR();
    double s = coord3d.getS();
    double t = coord3d.getT();
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 8 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(8));
    
    // d2H/dr2
    for (int i = 0; i < 8; i++) {
        deriv2[0][i] = 0.0;
    }
    
    // d2H/ds2
    for (int i = 0; i < 8; i++) {
        deriv2[1][i] = 0.0;
    }
    
    // d2H/dt2
    for (int i = 0; i < 8; i++) {
        deriv2[2][i] = 0.0;
    }
    
    // d2H/drds
    deriv2[3][0] =  0.125*(1 - t);
    deriv2[3][1] = -0.125*(1 - t);
    deriv2[3][2] =  0.125*(1 - t);
    deriv2[3][3] = -0.125*(1 - t);
    deriv2[3][4] =  0.125*(1 + t);
    deriv2[3][5] = -0.125*(1 + t);
    deriv2[3][6] =  0.125*(1 + t);
    deriv2[3][7] = -0.125*(1 + t);
    
    // d2H/drdt
    deriv2[4][0] =  0.125*(1 - s);
    deriv2[4][1] = -0.125*(1 - s);
    deriv2[4][2] = -0.125*(1 + s);
    deriv2[4][3] =  0.125*(1 + s);
    deriv2[4][4] = -0.125*(1 - s);
    deriv2[4][5] =  0.125*(1 - s);
    deriv2[4][6] =  0.125*(1 + s);
    deriv2[4][7] = -0.125*(1 + s);
    
    // d2H/dsdt
    deriv2[5][0] =  0.125*(1 - r);
    deriv2[5][1] =  0.125*(1 + r);
    deriv2[5][2] = -0.125*(1 + r);
    deriv2[5][3] = -0.125*(1 - r);
    deriv2[5][4] = -0.125*(1 - r);
    deriv2[5][5] = -0.125*(1 + r);
    deriv2[5][6] =  0.125*(1 + r);
    deriv2[5][7] =  0.125*(1 - r);
    
    return deriv2;
}