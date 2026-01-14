#include "RgShellElementTraits.h"
#include "RgElementTraits.h"
#include <assert.h>
#include <cmath>
#include "../RgGaussPoint.h"
#include "femcore/RgElementShapeStore.h"
#include "../ElementShape/RgSolidElementShape.h"
#include "elements/NaturalCoord.h"

using namespace RgFem;

//=============================================================================
// RgShellElementTraits
//=============================================================================

RgShellElementTraits::RgShellElementTraits(int ni, int ne, ElementShape es, ElementType et)
    : RgElementTraits(ni, ne, FE_ELEM_SHELL, es, et)
{
    // Initialize gaussPoints vector
    gaussPoints.resize(ni);
    
    Hr.resize(ni, ne);
    Hs.resize(ni, ne);
}

void RgShellElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);

    // calculate shape function values at gauss points
    for (int n = 0; n < m_nint; ++n)
    {
        NaturalCoord coord(gaussPoints[n].getR(), gaussPoints[n].getS());
        auto N = evalH(coord);
        for (int i = 0; i < m_neln; ++i)
            m_H(n, i) = N[i];
    }

    // calculate local derivatives of shape functions at gauss points
    for (int n = 0; n < m_nint; ++n)
    {
        NaturalCoord coord(gaussPoints[n].getR(), gaussPoints[n].getS());
        auto deriv = evalDeriv(coord);
        for (int i = 0; i < m_neln; ++i)
        {
            Hr(n, i) = deriv[0][i];  // dr derivatives
            Hs(n, i) = deriv[1][i];  // ds derivatives
        }
    }
}
std::vector<double> RgShellElementTraits::evalH(const NaturalCoord& coord)
{
    return m_shape->evalH(coord);
}

std::vector<std::vector<double>> RgShellElementTraits::evalDeriv(const NaturalCoord& coord)
{
    return m_shape->evalDeriv(coord);
}

std::vector<std::vector<double>> RgShellElementTraits::evalDeriv2(const NaturalCoord& coord)
{
    return m_shape->evalDeriv2(coord);
}

RgShellQuad4_::RgShellQuad4_(int ni, ElementType et) : RgShellElementTraits(ni, NELN, ET_QUAD4, et)
{
}

//=============================================================================
// RgShellQuad4G8
//=============================================================================

RgShellQuad4G8::RgShellQuad4G8() : RgShellQuad4_(NINT, FE_SHELL_QUAD4G8)
{
    // 4-point integration for in-plane, 2 for thickness
    // integration point coordinates
    const double a = 1.0 / sqrt(3.0);
    
    // Setup 8 integration points (4 in-plane x 2 through-thickness)
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            int idx = i * 2 + j;
            // in-plane coordinates
            gaussPoints[idx] = RgGaussPoint((i < 2) ? -a : a, (i % 2 == 0) ? -a : a, 0, 1.0); // simplified for shell
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellQuad4G8::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellQuad4G8::ni[NELN] = {};

//=============================================================================
// RgShellQuad4G12
//=============================================================================

RgShellQuad4G12::RgShellQuad4G12() : RgShellQuad4_(NINT, FE_SHELL_QUAD4G12)
{
    // 4-point integration for in-plane, 3 for thickness
    // integration point coordinates
    const double a = sqrt(0.6);
    
    // Setup 12 integration points (4 in-plane x 3 through-thickness)
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int idx = i * 3 + j;
            // in-plane coordinates
            gaussPoints[idx] = RgGaussPoint((i < 2) ? -a : a, (i % 2 == 0) ? -a : a, 0, 1.0); // simplified for shell
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellQuad4G12::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellQuad4G12::ni[NELN] = {};


RgShellTri3_::RgShellTri3_(int ni, ElementType et) : RgShellElementTraits(ni, NELN, ET_TRI3, et)
{
}

//=============================================================================
// RgShellTri3G6
//=============================================================================

RgShellTri3G6::RgShellTri3G6() : RgShellTri3_(NINT, FE_SHELL_TRI3G6)
{
    // 3-point integration for in-plane, 2 for thickness
    // Integration points for triangle with 3-point rule
    const double a = 1.0/6.0;
    const double b = 2.0/3.0;
    
    // Setup 6 integration points (3 in-plane x 2 through-thickness)
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            int idx = i * 2 + j;
            // in-plane coordinates
            if (i == 0) {
                gaussPoints[idx] = RgGaussPoint(a, a, 0, 0.5); // simplified for shell
            } else if (i == 1) {
                gaussPoints[idx] = RgGaussPoint(b, a, 0, 0.5);
            } else {
                gaussPoints[idx] = RgGaussPoint(a, b, 0, 0.5);
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellTri3G6::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellTri3G6::ni[NELN] = {};

//=============================================================================
// RgShellTri3G9
//=============================================================================

RgShellTri3G9::RgShellTri3G9() : RgShellTri3_(NINT, FE_SHELL_TRI3G9)
{
    // 3-point integration for in-plane, 3 for thickness
    // Integration points for triangle with 3-point rule
    const double a = 1.0/6.0;
    const double b = 2.0/3.0;
    
    // Setup 9 integration points (3 in-plane x 3 through-thickness)
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int idx = i * 3 + j;
            // in-plane coordinates
            if (i == 0) {
                gaussPoints[idx] = RgGaussPoint(a, a, 0, 0.5); // simplified for shell
            } else if (i == 1) {
                gaussPoints[idx] = RgGaussPoint(b, a, 0, 0.5);
            } else {
                gaussPoints[idx] = RgGaussPoint(a, b, 0, 0.5);
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellTri3G9::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellTri3G9::ni[NELN] = {};


RgShellQuad8_::RgShellQuad8_(int ni, ElementType et) : RgShellElementTraits(ni, NELN, ET_QUAD8, et)
{
}

//=============================================================================
// RgShellQuad8G18
//=============================================================================

RgShellQuad8G18::RgShellQuad8G18() : RgShellQuad8_(NINT, FE_SHELL_QUAD8G18)
{
    // 9-point integration for in-plane, 2 for thickness
    const double a = sqrt(0.6);
    const double w = 0.5555555555555556; // weight for gauss point
    const double w2 = 0.8888888888888888; // weight for center point
    
    // Setup 18 integration points (9 in-plane x 2 through-thickness)
    for (int i = 0; i < 9; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            int idx = i * 2 + j;
            // in-plane coordinates (using 3x3 gauss points)
            if (i == 0) { 
                gaussPoints[idx] = RgGaussPoint(-a, -a, 0, w * w); 
            } else if (i == 1) { 
                gaussPoints[idx] = RgGaussPoint(0, -a, 0, w2 * w); 
            } else if (i == 2) { 
                gaussPoints[idx] = RgGaussPoint(a, -a, 0, w * w); 
            } else if (i == 3) { 
                gaussPoints[idx] = RgGaussPoint(-a, 0, 0, w * w2); 
            } else if (i == 4) { 
                gaussPoints[idx] = RgGaussPoint(0, 0, 0, w2 * w2); 
            } else if (i == 5) { 
                gaussPoints[idx] = RgGaussPoint(a, 0, 0, w * w2); 
            } else if (i == 6) { 
                gaussPoints[idx] = RgGaussPoint(-a, a, 0, w * w); 
            } else if (i == 7) { 
                gaussPoints[idx] = RgGaussPoint(0, a, 0, w2 * w); 
            } else { 
                gaussPoints[idx] = RgGaussPoint(a, a, 0, w * w); 
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellQuad8G18::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellQuad8G18::ni[NELN] = {};

//=============================================================================
// RgShellQuad8G27
//=============================================================================

RgShellQuad8G27::RgShellQuad8G27() : RgShellQuad8_(NINT, FE_SHELL_QUAD8G27)
{
    // 9-point integration for in-plane, 3 for thickness
    const double a = sqrt(0.6);
    const double w = 0.5555555555555556; // weight for gauss point
    const double w2 = 0.8888888888888888; // weight for center point
    
    // Setup 27 integration points (9 in-plane x 3 through-thickness)
    for (int i = 0; i < 9; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int idx = i * 3 + j;
            // in-plane coordinates (using 3x3 gauss points)
            if (i == 0) { 
                gaussPoints[idx] = RgGaussPoint(-a, -a, 0, w * w); 
            } else if (i == 1) { 
                gaussPoints[idx] = RgGaussPoint(0, -a, 0, w2 * w); 
            } else if (i == 2) { 
                gaussPoints[idx] = RgGaussPoint(a, -a, 0, w * w); 
            } else if (i == 3) { 
                gaussPoints[idx] = RgGaussPoint(-a, 0, 0, w * w2); 
            } else if (i == 4) { 
                gaussPoints[idx] = RgGaussPoint(0, 0, 0, w2 * w2); 
            } else if (i == 5) { 
                gaussPoints[idx] = RgGaussPoint(a, 0, 0, w * w2); 
            } else if (i == 6) { 
                gaussPoints[idx] = RgGaussPoint(-a, a, 0, w * w); 
            } else if (i == 7) { 
                gaussPoints[idx] = RgGaussPoint(0, a, 0, w2 * w); 
            } else { 
                gaussPoints[idx] = RgGaussPoint(a, a, 0, w * w); 
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellQuad8G27::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellQuad8G27::ni[NELN] = {};


RgShellTri6_::RgShellTri6_(int ni, ElementType et) : RgShellElementTraits(ni, NELN, ET_TRI6, et)
{
}

//=============================================================================
// RgShellTri6G14
//=============================================================================

RgShellTri6G14::RgShellTri6G14() : RgShellTri6_(NINT, FE_SHELL_TRI6G14)
{
    // 7-point integration for in-plane, 2 for thickness
    // 7-point rule for triangle
    const double a = 0.059715871789770;
    const double b = 0.470142064105115;
    const double c = 0.797426985353087;
    const double d = 0.101286507323456;
    const double e = 0.333333333333333;
    const double w1 = 0.156892908754904;
    const double w2 = 0.125939180544827;
    const double w3 = 0.132394152788506;
    
    // Setup 14 integration points (7 in-plane x 2 through-thickness)
    for (int i = 0; i < 7; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            int idx = i * 2 + j;
            // in-plane coordinates
            if (i == 0) { 
                gaussPoints[idx] = RgGaussPoint(a, a, 0, w1); 
            } else if (i == 1) { 
                gaussPoints[idx] = RgGaussPoint(b, a, 0, w1); 
            } else if (i == 2) { 
                gaussPoints[idx] = RgGaussPoint(a, b, 0, w1); 
            } else if (i == 3) { 
                gaussPoints[idx] = RgGaussPoint(c, d, 0, w2); 
            } else if (i == 4) { 
                gaussPoints[idx] = RgGaussPoint(d, c, 0, w2); 
            } else if (i == 5) { 
                gaussPoints[idx] = RgGaussPoint(d, d, 0, w2); 
            } else { 
                gaussPoints[idx] = RgGaussPoint(e, e, 0, w3); 
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellTri6G14::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellTri6G14::ni[NELN] = {};

//=============================================================================
// RgShellTri6G21
//=============================================================================

RgShellTri6G21::RgShellTri6G21() : RgShellTri6_(NINT, FE_SHELL_TRI6G21)
{
    // 7-point integration for in-plane, 3 for thickness
    // 7-point rule for triangle
    const double a = 0.059715871789770;
    const double b = 0.470142064105115;
    const double c = 0.797426985353087;
    const double d = 0.101286507323456;
    const double e = 0.333333333333333;
    const double w1 = 0.156892908754904;
    const double w2 = 0.125939180544827;
    const double w3 = 0.132394152788506;
    
    // Setup 21 integration points (7 in-plane x 3 through-thickness)
    for (int i = 0; i < 7; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int idx = i * 3 + j;
            // in-plane coordinates
            if (i == 0) { 
                gaussPoints[idx] = RgGaussPoint(a, a, 0, w1); 
            } else if (i == 1) { 
                gaussPoints[idx] = RgGaussPoint(b, a, 0, w1); 
            } else if (i == 2) { 
                gaussPoints[idx] = RgGaussPoint(a, b, 0, w1); 
            } else if (i == 3) { 
                gaussPoints[idx] = RgGaussPoint(c, d, 0, w2); 
            } else if (i == 4) { 
                gaussPoints[idx] = RgGaussPoint(d, c, 0, w2); 
            } else if (i == 5) { 
                gaussPoints[idx] = RgGaussPoint(d, d, 0, w2); 
            } else { 
                gaussPoints[idx] = RgGaussPoint(e, e, 0, w3); 
            }
        }
    }
    
    init();
    
    // Compute inverse matrix for projection
    Matrix H(NINT, NELN);
    for (int i = 0; i < NINT; ++i)
    {
        double r = gaussPoints[i].getR();
        double s = gaussPoints[i].getS();
        auto h = evalH(NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Hi = H.inverse();
}

void RgShellTri6G21::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Hi[j][k] * ai[k];
        }
    }
}

int RgShellTri6G21::ni[NELN] = {};