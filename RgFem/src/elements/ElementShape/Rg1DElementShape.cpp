#include "Rg1DElementShape.h"



// Implementation for 2-node linear 1D element
std::vector<double> Rg1D2NodeElementShape::evalH(const NaturalCoord& coord)
{
    double r = coord.getR();
    std::vector<double> H(2);
    
    // Linear 1D shape functions in natural coordinates [-1, 1]
    H[0] = 0.5 * (1.0 - r);  // Node 1 at r = -1
    H[1] = 0.5 * (1.0 + r);  // Node 2 at r = +1
    
    return H;
}

std::vector<std::vector<double>> Rg1D2NodeElementShape::evalDeriv(const NaturalCoord& coord)
{
    // For 1D element, derivatives are with respect to the single natural coordinate
    std::vector<std::vector<double>> deriv(1, std::vector<double>(2));
    
    // Derivatives of shape functions with respect to r
    deriv[0][0] = -0.5;  // dN1/dr
    deriv[0][1] = 0.5;   // dN2/dr
    
    return deriv;
}

std::vector<std::vector<double>> Rg1D2NodeElementShape::evalDeriv2(const NaturalCoord& coord)
{
    // For 1D element, second derivatives
    std::vector<std::vector<double>> deriv2(1, std::vector<double>(2));
    
    // Second derivatives of shape functions with respect to r
    deriv2[0][0] = 0.0;  // d2N1/dr2
    deriv2[0][1] = 0.0;  // d2N2/dr2
    
    return deriv2;
}

// Implementation for 3-node quadratic 1D element
std::vector<double> Rg1D3NodeElementShape::evalH(const NaturalCoord& coord)
{
    double r = coord.getR();
    std::vector<double> H(3);
    
    // Quadratic 1D shape functions in natural coordinates [-1, 1]
    // Node 1 at r = -1, Node 2 at r = 0, Node 3 at r = +1
    H[0] = 0.5 * r * (r - 1.0);   // N1
    H[1] = (1.0 - r) * (1.0 + r); // N2
    H[2] = 0.5 * r * (r + 1.0);   // N3
    
    return H;
}

std::vector<std::vector<double>> Rg1D3NodeElementShape::evalDeriv(const NaturalCoord& coord)
{
    // For 1D element, derivatives are with respect to the single natural coordinate
    std::vector<std::vector<double>> deriv(1, std::vector<double>(3));
    
    double r = coord.getR();
    
    // Derivatives of shape functions with respect to r
    deriv[0][0] = r - 0.5;      // dN1/dr
    deriv[0][1] = -2.0 * r;     // dN2/dr
    deriv[0][2] = r + 0.5;      // dN3/dr
    
    return deriv;
}

std::vector<std::vector<double>> Rg1D3NodeElementShape::evalDeriv2(const NaturalCoord& coord)
{
    // For 1D element, second derivatives
    std::vector<std::vector<double>> deriv2(1, std::vector<double>(3));
    
    // Second derivatives of shape functions with respect to r
    deriv2[0][0] = 1.0;   // d2N1/dr2
    deriv2[0][1] = -2.0;  // d2N2/dr2
    deriv2[0][2] = 1.0;   // d2N3/dr2
    
    return deriv2;
}

// Implementation for 4-node cubic 1D element
std::vector<double> Rg1D4NodeElementShape::evalH(const NaturalCoord& coord)
{
    double r = coord.getR();
    std::vector<double> H(4);
    
    // Cubic 1D shape functions in natural coordinates [-1, 1]
    // Node 1 at r = -1, Node 2 at -1/3, Node 3 at 1/3, Node 4 at r = +1
    H[0] = -0.5 * (r - 0.0) * (r - 1.0/3.0) * (r - 1.0) / ((-1.0/3.0) * (2.0/3.0) * 2.0); // N1 at r=-1
    H[1] = (r + 1.0) * (r - 1.0/3.0) * (r - 1.0) / ((4.0/3.0) * (-2.0/3.0) * (-4.0/3.0)); // N2 at r=-1/3
    H[2] = (r + 1.0) * (r + 1.0/3.0) * (r - 1.0) / ((2.0/3.0) * (4.0/3.0) * (-2.0/3.0)); // N3 at r=1/3
    H[3] = (r + 1.0) * (r + 1.0/3.0) * (r - 0.0) / (2.0 * (4.0/3.0) * (2.0/3.0)); // N4 at r=1
    
    // Simplified Lagrange cubic shape functions
    // N1 = -0.5 * r * (r-1/3) * (r-1) / (2/9)
    // N2 = (r+1) * (r-1/3) * (r-1) / (64/27)
    // N3 = (r+1) * (r+1/3) * (r-1) / (-64/27)
    // N4 = 0.5 * (r+1) * (r+1/3) * r / (2/9)
    
    // More direct cubic Lagrange shape functions:
    // At r = [-1, -1/3, 1/3, 1]
    double r1 = -1.0, r2 = -1.0/3.0, r3 = 1.0/3.0, r4 = 1.0;
    
    // N1(r) = [(r-r2)(r-r3)(r-r4)] / [(r1-r2)(r1-r3)(r1-r4)]
    double den1 = (r1-r2)*(r1-r3)*(r1-r4);
    H[0] = ((r-r2)*(r-r3)*(r-r4)) / den1;
    
    // N2(r) = [(r-r1)(r-r3)(r-r4)] / [(r2-r1)(r2-r3)(r2-r4)]
    double den2 = (r2-r1)*(r2-r3)*(r2-r4);
    H[1] = ((r-r1)*(r-r3)*(r-r4)) / den2;
    
    // N3(r) = [(r-r1)(r-r2)(r-r4)] / [(r3-r1)(r3-r2)(r3-r4)]
    double den3 = (r3-r1)*(r3-r2)*(r3-r4);
    H[2] = ((r-r1)*(r-r2)*(r-r4)) / den3;
    
    // N4(r) = [(r-r1)(r-r2)(r-r3)] / [(r4-r1)(r4-r2)(r4-r3)]
    double den4 = (r4-r1)*(r4-r2)*(r4-r3);
    H[3] = ((r-r1)*(r-r2)*(r-r3)) / den4;
    
    return H;
}

std::vector<std::vector<double>> Rg1D4NodeElementShape::evalDeriv(const NaturalCoord& coord)
{
    // For 1D element, derivatives are with respect to the single natural coordinate
    std::vector<std::vector<double>> deriv(1, std::vector<double>(4));
    
    double r = coord.getR();
    
    // Derivatives of cubic shape functions
    double r1 = -1.0, r2 = -1.0/3.0, r3 = 1.0/3.0, r4 = 1.0;
    double den1 = (r1-r2)*(r1-r3)*(r1-r4);
    double den2 = (r2-r1)*(r2-r3)*(r2-r4);
    double den3 = (r3-r1)*(r3-r2)*(r3-r4);
    double den4 = (r4-r1)*(r4-r2)*(r4-r3);
    
    // dN1/dr = [ (r-r3)(r-r4) + (r-r2)(r-r4) + (r-r2)(r-r3) ] / den1
    deriv[0][0] = ((r-r3)*(r-r4) + (r-r2)*(r-r4) + (r-r2)*(r-r3)) / den1;
    
    // dN2/dr = [ (r-r3)(r-r4) + (r-r1)(r-r4) + (r-r1)(r-r3) ] / den2
    deriv[0][1] = ((r-r3)*(r-r4) + (r-r1)*(r-r4) + (r-r1)*(r-r3)) / den2;
    
    // dN3/dr = [ (r-r2)(r-r4) + (r-r1)(r-r4) + (r-r1)(r-r2) ] / den3
    deriv[0][2] = ((r-r2)*(r-r4) + (r-r1)*(r-r4) + (r-r1)*(r-r2)) / den3;
    
    // dN4/dr = [ (r-r2)(r-r3) + (r-r1)(r-r3) + (r-r1)(r-r2) ] / den4
    deriv[0][3] = ((r-r2)*(r-r3) + (r-r1)*(r-r3) + (r-r1)*(r-r2)) / den4;
    
    return deriv;
}

std::vector<std::vector<double>> Rg1DHermiteElementShape::evalDeriv(const NaturalCoord& coord)
{
    // For 1D element, derivatives are with respect to the single natural coordinate
    std::vector<std::vector<double>> deriv(1, std::vector<double>(4));
    
    double r = coord.getR();
    
    // First derivatives of Hermite shape functions
    deriv[0][0] = 0.25 * (-3.0 + 3.0*r*r);     // dN1/dr
    deriv[0][1] = 0.25 * (3.0 - 3.0*r*r);      // dN2/dr
    deriv[0][2] = 0.125 * (-1.0 - 2.0*r - 3.0*r*r);  // dN3/dr
    deriv[0][3] = 0.125 * (1.0 - 2.0*r + 3.0*r*r);   // dN4/dr
    
    return deriv;
}

std::vector<std::vector<double>> Rg1D4NodeElementShape::evalDeriv2(const NaturalCoord& coord)
{
    // For 1D element, second derivatives
    std::vector<std::vector<double>> deriv2(1, std::vector<double>(4));
    
    double r = coord.getR();
    
    // Second derivatives of cubic shape functions
    double r1 = -1.0, r2 = -1.0/3.0, r3 = 1.0/3.0, r4 = 1.0;
    double den1 = (r1-r2)*(r1-r3)*(r1-r4);
    double den2 = (r2-r1)*(r2-r3)*(r2-r4);
    double den3 = (r3-r1)*(r3-r2)*(r3-r4);
    double den4 = (r4-r1)*(r4-r2)*(r4-r3);
    
    // d2N1/dr2 = [2(r-r3) + 2(r-r2) + 2(r-r4)] / den1
    deriv2[0][0] = (2.0*(r-r3) + 2.0*(r-r2) + 2.0*(r-r4)) / den1;
    
    // d2N2/dr2 = [2(r-r3) + 2(r-r1) + 2(r-r4)] / den2
    deriv2[0][1] = (2.0*(r-r3) + 2.0*(r-r1) + 2.0*(r-r4)) / den2;
    
    // d2N3/dr2 = [2(r-r2) + 2(r-r1) + 2(r-r4)] / den3
    deriv2[0][2] = (2.0*(r-r2) + 2.0*(r-r1) + 2.0*(r-r4)) / den3;
    
    // d2N4/dr2 = [2(r-r2) + 2(r-r1) + 2(r-r3)] / den4
    deriv2[0][3] = (2.0*(r-r2) + 2.0*(r-r1) + 2.0*(r-r3)) / den4;
    
    return deriv2;
}

std::vector<std::vector<double>> Rg1DHermiteElementShape::evalDeriv2(const NaturalCoord& coord)
{
    // For 1D element, second derivatives
    std::vector<std::vector<double>> deriv2(1, std::vector<double>(4));
    
    double r = coord.getR();
    
    // Second derivatives of Hermite shape functions
    deriv2[0][0] = 0.25 * (6.0*r);              // d2N1/dr2
    deriv2[0][1] = 0.25 * (-6.0*r);             // d2N2/dr2
    deriv2[0][2] = 0.125 * (-2.0 - 6.0*r);      // d2N3/dr2
    deriv2[0][3] = 0.125 * (-2.0 + 6.0*r);      // d2N4/dr2
    
    return deriv2;
}

