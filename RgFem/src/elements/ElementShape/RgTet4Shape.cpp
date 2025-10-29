#include "elements/ElementShape/RgTet4Shape.h"
#include <vector>
#include <stdexcept>

std::vector<double> RgTet4Shape::evalH(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgTet4Shape requires NaturalCoord3d coordinates");
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

std::vector<std::vector<double>> RgTet4Shape::evalDeriv(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgTet4Shape requires NaturalCoord3d coordinates");
    }
    
    // Return derivatives in the format [dH/dr, dH/ds, dH/dt]
    // Each derivative is a vector of size 4 (one for each node)
    std::vector<std::vector<double>> deriv(3, std::vector<double>(4));
    
    // dH/dr
    deriv[0][0] = -1; deriv[0][1] = 1; deriv[0][2] = 0; deriv[0][3] = 0;
    
    // dH/ds
    deriv[1][0] = -1; deriv[1][1] = 0; deriv[1][2] = 1; deriv[1][3] = 0;
    
    // dH/dt
    deriv[2][0] = -1; deriv[2][1] = 0; deriv[2][2] = 0; deriv[2][3] = 1;
    
    return deriv;
}

std::vector<std::vector<double>> RgTet4Shape::evalDeriv2(const NaturalCoord& coord)
{
    // Check if the coordinate is of the correct type
    const NaturalCoord3d* coord3d = dynamic_cast<const NaturalCoord3d*>(&coord);
    if (!coord3d) {
        throw std::invalid_argument("RgTet4Shape requires NaturalCoord3d coordinates");
    }
    
    // Return second derivatives in the format:
    // [d2H/dr2, d2H/ds2, d2H/dt2, d2H/drds, d2H/dsdt, d2H/drdt]
    // Each derivative is a vector of size 4 (one for each node)
    std::vector<std::vector<double>> deriv2(6, std::vector<double>(4, 0.0));
    
    // All second derivatives are zero for linear shape functions
    // So we just return the initialized zero vectors
    
    return deriv2;
}