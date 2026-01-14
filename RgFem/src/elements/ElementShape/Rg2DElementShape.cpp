#include "Rg2DElementShape.h"
#include "elements/NaturalCoord.h"
#include <vector>



//=============================================================================
//              Q U A D 4
//=============================================================================

//-----------------------------------------------------------------------------
std::vector<double> FEQuad4::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(4);
    double r = coord.getR();
    double s = coord.getS();
    
    H[0] = 0.25*(1 - r)*(1 - s);
    H[1] = 0.25*(1 + r)*(1 - s);
    H[2] = 0.25*(1 + r)*(1 + s);
    H[3] = 0.25*(1 - r)*(1 + s);
    
    return H;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FEQuad4::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 4;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] = -0.25*(1 - s); result[1][0] = -0.25*(1 - r);
    result[0][1] =  0.25*(1 - s); result[1][1] = -0.25*(1 + r);
    result[0][2] =  0.25*(1 + s); result[1][2] =  0.25*(1 + r);
    result[0][3] = -0.25*(1 + s); result[1][3] =  0.25*(1 - r);
    
    return result;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FEQuad4::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 4;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] = 0; result[2][0] =  0.25; result[1][0] = 0;
    result[0][1] = 0; result[2][1] = -0.25; result[1][1] = 0;
    result[0][2] = 0; result[2][2] =  0.25; result[1][2] = 0;
    result[0][3] = 0; result[2][3] = -0.25; result[1][3] = 0;
    
    return result;
}

//=============================================================================
//              Q U A D 8
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
std::vector<double> FEQuad8::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(8);
    double r = coord.getR();
    double s = coord.getS();
    
    H[4] = 0.5*(1 - r*r)*(1 - s);
    H[5] = 0.5*(1 - s*s)*(1 + r);
    H[6] = 0.5*(1 - r*r)*(1 + s);
    H[7] = 0.5*(1 - s*s)*(1 - r);

    H[0] = 0.25*(1 - r)*(1 - s) - 0.5*(H[4] + H[7]);
    H[1] = 0.25*(1 + r)*(1 - s) - 0.5*(H[4] + H[5]);
    H[2] = 0.25*(1 + r)*(1 + s) - 0.5*(H[5] + H[6]);
    H[3] = 0.25*(1 - r)*(1 + s) - 0.5*(H[6] + H[7]);
    
    return H;
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
std::vector<std::vector<double>> FEQuad8::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 8;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][4] = -r*(1 - s);
    result[0][5] = 0.5*(1 - s*s);
    result[0][6] = -r*(1 + s);
    result[0][7] = -0.5*(1 - s*s);

    result[0][0] = -0.25*(1 - s) - 0.5*(result[0][4] + result[0][7]);
    result[0][1] = 0.25*(1 - s) - 0.5*(result[0][4] + result[0][5]);
    result[0][2] = 0.25*(1 + s) - 0.5*(result[0][5] + result[0][6]);
    result[0][3] = -0.25*(1 + s) - 0.5*(result[0][6] + result[0][7]);

    result[1][4] = -0.5*(1 - r*r);
    result[1][5] = -s*(1 + r);
    result[1][6] = 0.5*(1 - r*r);
    result[1][7] = -s*(1 - r);

    result[1][0] = -0.25*(1 - r) - 0.5*(result[1][4] + result[1][7]);
    result[1][1] = -0.25*(1 + r) - 0.5*(result[1][4] + result[1][5]);
    result[1][2] = 0.25*(1 + r) - 0.5*(result[1][5] + result[1][6]);
    result[1][3] = 0.25*(1 - r) - 0.5*(result[1][6] + result[1][7]);
    
    return result;
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
std::vector<std::vector<double>> FEQuad8::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 8;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][4] = -(1 - s);
    result[0][5] = 0.0;
    result[0][6] = -(1 + s);
    result[0][7] = 0.0;

    result[2][4] = r;
    result[2][5] = -s;
    result[2][6] = -r;
    result[2][7] = s;

    result[1][4] = 0.0;
    result[1][5] = -(1 + r);
    result[1][6] = 0.0;
    result[1][7] = -(1 - r);

    result[0][0] = -0.5*(result[0][4] + result[0][7]);
    result[0][1] = -0.5*(result[0][4] + result[0][5]);
    result[0][2] = -0.5*(result[0][5] + result[0][6]);
    result[0][3] = -0.5*(result[0][6] + result[0][7]);

    result[2][0] = 0.25 - 0.5*(result[2][4] + result[2][7]);
    result[2][1] = -0.25 - 0.5*(result[2][4] + result[2][5]);
    result[2][2] = 0.25 - 0.5*(result[2][5] + result[2][6]);
    result[2][3] = -0.25 - 0.5*(result[2][6] + result[2][7]);

    result[1][0] = -0.5*(result[1][4] + result[1][7]);
    result[1][1] = -0.5*(result[1][4] + result[1][5]);
    result[1][2] = -0.5*(result[1][5] + result[1][6]);
    result[1][3] = -0.5*(result[1][6] + result[1][7]);
    
    return result;
}

//=============================================================================
//              Q U A D 9
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
std::vector<double> FEQuad9::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(9);
    double r = coord.getR();
    double s = coord.getS();
    
    double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
    double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };

    H[0] = R[0] * S[0];
    H[1] = R[1] * S[0];
    H[2] = R[1] * S[1];
    H[3] = R[0] * S[1];
    H[4] = R[2] * S[0];
    H[5] = R[1] * S[2];
    H[6] = R[2] * S[1];
    H[7] = R[0] * S[2];
    H[8] = R[2] * S[2];
    
    return H;
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
std::vector<std::vector<double>> FEQuad9::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 9;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
    double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
    double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
    double DS[3] = { s - 0.5, s + 0.5, -2.0*s };

    result[0][0] = DR[0] * S[0];
    result[0][1] = DR[1] * S[0];
    result[0][2] = DR[1] * S[1];
    result[0][3] = DR[0] * S[1];
    result[0][4] = DR[2] * S[0];
    result[0][5] = DR[1] * S[2];
    result[0][6] = DR[2] * S[1];
    result[0][7] = DR[0] * S[2];
    result[0][8] = DR[2] * S[2];

    result[1][0] = R[0] * DS[0];
    result[1][1] = R[1] * DS[0];
    result[1][2] = R[1] * DS[1];
    result[1][3] = R[0] * DS[1];
    result[1][4] = R[2] * DS[0];
    result[1][5] = R[1] * DS[2];
    result[1][6] = R[2] * DS[1];
    result[1][7] = R[0] * DS[2];
    result[1][8] = R[2] * DS[2];
    
    return result;
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
std::vector<std::vector<double>> FEQuad9::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 9;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
    double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
    double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
    double DS[3] = { s - 0.5, s + 0.5, -2.0*s };
    double DDR[3] = { 1.0, 1.0, -2.0 };
    double DDS[3] = { 1.0, 1.0, -2.0 };

    result[0][0] = DDR[0] * S[0]; result[2][0] = DR[0] * DS[0]; result[1][0] = R[0] * DDS[0];
    result[0][1] = DDR[1] * S[0]; result[2][1] = DR[1] * DS[0]; result[1][1] = R[1] * DDS[0];
    result[0][2] = DDR[1] * S[1]; result[2][2] = DR[1] * DS[1]; result[1][2] = R[1] * DDS[1];
    result[0][3] = DDR[0] * S[1]; result[2][3] = DR[0] * DS[1]; result[1][3] = R[0] * DDS[1];
    result[0][4] = DDR[2] * S[0]; result[2][4] = DR[2] * DS[0]; result[1][4] = R[2] * DDS[0];
    result[0][5] = DDR[1] * S[2]; result[2][5] = DR[1] * DS[2]; result[1][5] = R[1] * DDS[2];
    result[0][6] = DDR[2] * S[1]; result[2][6] = DR[2] * DS[1]; result[1][6] = R[2] * DDS[1];
    result[0][7] = DDR[0] * S[2]; result[2][7] = DR[0] * DS[2]; result[1][7] = R[0] * DDS[2];
    result[0][8] = DDR[2] * S[2]; result[2][8] = DR[2] * DS[2]; result[1][8] = R[2] * DDS[2];
    
    return result;
}

//=============================================================================
//              T R I 3
//=============================================================================

//-----------------------------------------------------------------------------
std::vector<double> FETri3::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(3);
    double r = coord.getR();
    double s = coord.getS();
    
    H[0] = 1.0 - r - s;
    H[1] = r;
    H[2] = s;
    
    return H;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri3::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 3;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] = -1; result[1][0] = -1;
    result[0][1] =  1; result[1][1] =  0;
    result[0][2] =  0; result[1][2] =  1;
    
    return result;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri3::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 3;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] = 0; result[2][0] = 0; result[1][0] = 0;
    result[0][1] = 0; result[2][1] = 0; result[1][1] = 0;
    result[0][2] = 0; result[2][2] = 0; result[1][2] = 0;
    
    return result;
}

//=============================================================================
//              T R I 6
//=============================================================================

//-----------------------------------------------------------------------------
std::vector<double> FETri6::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(6);
    double r = coord.getR();
    double s = coord.getS();
    
    double r1 = 1.0 - r - s;
    double r2 = r;
    double r3 = s;

    H[0] = r1*(2.0*r1 - 1.0);
    H[1] = r2*(2.0*r2 - 1.0);
    H[2] = r3*(2.0*r3 - 1.0);
    H[3] = 4.0*r1*r2;
    H[4] = 4.0*r2*r3;
    H[5] = 4.0*r3*r1;
    
    return H;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri6::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 6;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] = -3.0 + 4.0*r + 4.0*s;
    result[0][1] = 4.0*r - 1.0;
    result[0][2] = 0.0;
    result[0][3] = 4.0 - 8.0*r - 4.0*s;
    result[0][4] = 4.0*s;
    result[0][5] = -4.0*s;

    result[1][0] = -3.0 + 4.0*s + 4.0*r;
    result[1][1] = 0.0;
    result[1][2] = 4.0*s - 1.0;
    result[1][3] = -4.0*r;
    result[1][4] = 4.0*r;
    result[1][5] = 4.0 - 8.0*s - 4.0*r;
    
    return result;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri6::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 6;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    result[0][0] =  4.0; result[2][0] =  4.0; result[1][0] =  4.0;
    result[0][1] =  4.0; result[2][1] =  0.0; result[1][1] =  0.0;
    result[0][2] =  0.0; result[2][2] =  0.0; result[1][2] =  4.0;
    result[0][3] = -8.0; result[2][3] = -4.0; result[1][3] =  0.0;
    result[0][4] =  0.0; result[2][4] =  4.0; result[1][4] =  0.0;
    result[0][5] =  0.0; result[2][5] = -4.0; result[1][5] = -8.0;
    
    return result;
}

//=============================================================================
//              T R I 7
//=============================================================================

//-----------------------------------------------------------------------------
std::vector<double> FETri7::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(7);
    double r = coord.getR();
    double s = coord.getS();
    
    double r1 = 1.0 - r - s;
    double r2 = r;
    double r3 = s;

    H[6] = 27.0*r1*r2*r3;
    H[0] = r1*(2.0*r1 - 1.0) + H[6] / 9.0;
    H[1] = r2*(2.0*r2 - 1.0) + H[6] / 9.0;
    H[2] = r3*(2.0*r3 - 1.0) + H[6] / 9.0;
    H[3] = 4.0*r1*r2 - 4.0*H[6] / 9.0;
    H[4] = 4.0*r2*r3 - 4.0*H[6] / 9.0;
    H[5] = 4.0*r3*r1 - 4.0*H[6] / 9.0;
    
    return H;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri7::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 7;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    double Hr[7];
    double Hs[7];
    
    Hr[6] = 27.0*s*(1.0 - 2.0*r - s);
    Hr[0] = -3.0 + 4.0*r + 4.0*s + Hr[6] / 9.0;
    Hr[1] = 4.0*r - 1.0 + Hr[6] / 9.0;
    Hr[2] = 0.0 + Hr[6] / 9.0;
    Hr[3] = 4.0 - 8.0*r - 4.0*s - 4.0*Hr[6] / 9.0;
    Hr[4] = 4.0*s - 4.0*Hr[6] / 9.0;
    Hr[5] = -4.0*s - 4.0*Hr[6] / 9.0;

    Hs[6] = 27.0*r*(1.0 - r - 2.0*s);
    Hs[0] = -3.0 + 4.0*s + 4.0*r + Hs[6] / 9.0;
    Hs[1] = 0.0 + Hs[6] / 9.0;
    Hs[2] = 4.0*s - 1.0 + Hs[6] / 9.0;
    Hs[3] = -4.0*r - 4.0*Hs[6] / 9.0;
    Hs[4] = 4.0*r - 4.0*Hs[6] / 9.0;
    Hs[5] = 4.0 - 8.0*s - 4.0*r - 4.0*Hs[6] / 9.0;
    
    result[0][0] = Hr[0]; result[1][0] = Hs[0];
    result[0][1] = Hr[1]; result[1][1] = Hs[1];
    result[0][2] = Hr[2]; result[1][2] = Hs[2];
    result[0][3] = Hr[3]; result[1][3] = Hs[3];
    result[0][4] = Hr[4]; result[1][4] = Hs[4];
    result[0][5] = Hr[5]; result[1][5] = Hs[5];
    result[0][6] = Hr[6]; result[1][6] = Hs[6];
    
    return result;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri7::evalDeriv2(const NaturalCoord& coord)
{
    int nNodes = 7;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    double r = coord.getR();
    double s = coord.getS();
    
    double Hrr[7];
    double Hss[7];
    double Hrs[7];
    
    Hrr[6] = -54.0*s;
    Hss[6] = -54.0*r;
    Hrs[6] = 27.0*(1.0 - 2.0*r - 2.0*s);

    Hrr[0] = 4.0 + Hrr[6] / 9.0; Hrs[0] = 4.0 + Hrs[6] / 9.0; Hss[0] = 4.0 + Hss[6] / 9.0;
    Hrr[1] = 4.0 + Hrr[6] / 9.0; Hrs[1] = 0.0 + Hrs[6] / 9.0; Hss[1] = 0.0 + Hss[6] / 9.0;
    Hrr[2] = 0.0 + Hrr[6] / 9.0; Hrs[2] = 0.0 + Hrs[6] / 9.0; Hss[2] = 4.0 + Hss[6] / 9.0;
    Hrr[3] = -8.0 - 4.0*Hrr[6] / 9.0; Hrs[3] = -4.0 - 4.0*Hrs[6] / 9.0; Hss[3] = 0.0 - 4.0*Hss[6] / 9.0;
    Hrr[4] = 0.0 - 4.0*Hrr[6] / 9.0; Hrs[4] = 4.0 - 4.0*Hrs[6] / 9.0; Hss[4] = 0.0 - 4.0*Hss[6] / 9.0;
    Hrr[5] = 0.0 - 4.0*Hrr[6] / 9.0; Hrs[5] = -4.0 - 4.0*Hrs[6] / 9.0; Hss[5] = -8.0 - 4.0*Hss[6] / 9.0;
    
    result[0][0] = Hrr[0]; result[2][0] = Hrs[0]; result[1][0] = Hss[0];
    result[0][1] = Hrr[1]; result[2][1] = Hrs[1]; result[1][1] = Hss[1];
    result[0][2] = Hrr[2]; result[2][2] = Hrs[2]; result[1][2] = Hss[2];
    result[0][3] = Hrr[3]; result[2][3] = Hrs[3]; result[1][3] = Hss[3];
    result[0][4] = Hrr[4]; result[2][4] = Hrs[4]; result[1][4] = Hss[4];
    result[0][5] = Hrr[5]; result[2][5] = Hrs[5]; result[1][5] = Hss[5];
    result[0][6] = Hrr[6]; result[2][6] = Hrs[6]; result[1][6] = Hss[6];
    
    return result;
}

//=============================================================================
//              T R I 10
//=============================================================================

//-----------------------------------------------------------------------------
std::vector<double> FETri10::evalH(const NaturalCoord& coord)
{
    std::vector<double> H(10);
    double r = coord.getR();
    double s = coord.getS();
    
    double L1 = 1.0 - r - s;
    double L2 = r;
    double L3 = s;

    H[0] = 0.5*(3 * L1 - 1)*(3 * L1 - 2)*L1;
    H[1] = 0.5*(3 * L2 - 1)*(3 * L2 - 2)*L2;
    H[2] = 0.5*(3 * L3 - 1)*(3 * L3 - 2)*L3;
    H[3] = 4.5*(3 * L1 - 1)*L1*L2;
    H[4] = 4.5*(3 * L2 - 1)*L1*L2;
    H[5] = 4.5*(3 * L2 - 1)*L2*L3;
    H[6] = 4.5*(3 * L3 - 1)*L2*L3;
    H[7] = 4.5*(3 * L1 - 1)*L1*L3;
    H[8] = 4.5*(3 * L3 - 1)*L1*L3;
    H[9] = 27.*L1*L2*L3;
    
    return H;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri10::evalDeriv(const NaturalCoord& coord)
{
    int nNodes = 10;
    std::vector<std::vector<double>> result(2, std::vector<double>(nNodes)); // 2 for d/dr, d/ds
    double r = coord.getR();
    double s = coord.getS();
    
    double L1 = 1.0 - r - s;
    double L2 = r;
    double L3 = s;

    result[0][0] = -3. / 2.*(3 * L1 - 2)*L1 - 3. / 2.*(3 * L1 - 1)*L1 - 0.5*(3 * L1 - 1)*(3 * L1 - 2);
    result[0][1] = 3. / 2.*(3 * L2 - 2)*L2 + 3. / 2.*(3 * L2 - 1)*L2 + 0.5*(3 * L2 - 1)*(3 * L2 - 2);
    result[0][2] = 0.0;
    result[0][3] = -27. / 2.*L1*L2 - 9. / 2.*(3 * L1 - 1)*L2 + 9. / 2.*(3 * L1 - 1)*L1;
    result[0][4] = 27. / 2.*L1*L2 - 9. / 2.*(3 * L2 - 1)*L2 + 9. / 2.*(3 * L2 - 1)*L1;
    result[0][5] = 27. / 2.*L2*L3 + 9. / 2.*(3 * L2 - 1)*L3;
    result[0][6] = 9. / 2.*(3 * L3 - 1)*L3;
    result[0][7] = -27. / 2.*L1*L3 - 9. / 2.*(3 * L1 - 1)*L3;
    result[0][8] = -9. / 2.*(3 * L3 - 1)*L3;
    result[0][9] = -27.*L2*L3 + 27.*L1*L3;

    result[1][0] = -3. / 2.*(3 * L1 - 2)*L1 - 3. / 2.*(3 * L1 - 1)*L1 - 0.5*(3 * L1 - 1)*(3 * L1 - 2);
    result[1][1] = 0.0;
    result[1][2] = 3. / 2.*(3 * L3 - 2)*L3 + 3. / 2.*(3 * L3 - 1)*L3 + 0.5*(3 * L3 - 1)*(3 * L3 - 2);
    result[1][3] = -27. / 2.*L1*L2 - 9. / 2.*(3 * L1 - 1)*L2;
    result[1][4] = -9. / 2.*(3 * L2 - 1)*L2;
    result[1][5] = 9. / 2.*(3 * L2 - 1)*L2;
    result[1][6] = 27. / 2.*L2*L3 + 9. / 2.*(3 * L3 - 1)*L2;
    result[1][7] = -27. / 2.*L1*L3 - 9. / 2.*(3 * L1 - 1)*L3 + 9. / 2.*(3 * L1 - 1)*L1;
    result[1][8] = 27. / 2.*L1*L3 - 9. / 2.*(3 * L3 - 1)*L3 + 9. / 2.*(3 * L3 - 1)*L1;
    result[1][9] = -27.*L2*L3 + 27.*L1*L2;
    
    return result;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<double>> FETri10::evalDeriv2(const NaturalCoord& coord)
{
    // TODO: Implement this
    int nNodes = 10;
    std::vector<std::vector<double>> result(3, std::vector<double>(nNodes)); // 3 for d2/dr2, d2/ds2, d2/drdt
    return result;
}