#include "Rg2DElementTraits.h"
#include "RgElementTraits.h"
#include <assert.h>
#include <cmath>
#include "../RgGaussPoint.h"
#include "femcore/RgElementShapeStore.h"
#include "../ElementShape/RgSolidElementShape.h"
#include "elements/NaturalCoord.h"

using namespace RgFem;

// Rg2DElementTraits implementation
Rg2DElementTraits::Rg2DElementTraits(int ni, int ne, ElementShape es, ElementType et)
    : RgElementTraits(ni, ne, FE_ELEM_SOLID_2D, es, et)
{
    gaussPoints.resize(ni);

    Gr.resize(ni, ne);
    Gs.resize(ni, ne);
    
    Grr.resize(ni, ne);
    Gsr.resize(ni, ne);
    
    Grs.resize(ni, ne);
    Gss.resize(ni, ne);
}

//-----------------------------------------------------------------------------
//! Initialize the 2D element traits data variables.
//
void Rg2DElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);

    // evaluate shape functions
    std::vector<double> N(m_neln);
    for (int n=0; n<m_nint; ++n)
    {
        shape(N.data(), gaussPoints[n].getR(), gaussPoints[n].getS());
        for (int i=0; i<m_neln; ++i) m_H(n, i) = N[i];
    }
    
    // evaluate shape function derivatives
    std::vector<double> Nr(m_neln);
    std::vector<double> Ns(m_neln);
    for (int n=0; n<m_nint; ++n)
    {
        shape_deriv(Nr.data(), Ns.data(), gaussPoints[n].getR(), gaussPoints[n].getS());
        for (int i=0; i<m_neln; ++i)
        {
            Gr(n, i) = Nr[i];
            Gs(n, i) = Ns[i];
        }
    }
}

//-----------------------------------------------------------------------------
// Implementation of modern interface using legacy interface
//
std::vector<double> Rg2DElementTraits::evalH(const RgFem::NaturalCoord& coord)
{
    const int NELN = m_neln;
    std::vector<double> N(NELN);
    shape(N.data(), coord.getR(), coord.getS());
    
    std::vector<double> result(NELN);
    for (int i = 0; i < NELN; ++i) {
        result[i] = N[i];
    }
    
    return result;
}

std::vector<std::vector<double>> Rg2DElementTraits::evalDeriv(const RgFem::NaturalCoord& coord)
{
    const int NELN = m_neln;
    std::vector<double> Nr(NELN), Ns(NELN);
    shape_deriv(Nr.data(), Ns.data(), coord.getR(), coord.getS());
    
    std::vector<std::vector<double>> result(2, std::vector<double>(NELN));
    for (int i = 0; i < NELN; ++i) {
        result[0][i] = Nr[i];  // dr derivatives
        result[1][i] = Ns[i];  // ds derivatives
    }
    
    return result;
}

//-----------------------------------------------------------------------------

//=============================================================================
//                          Rg2DTri3 
//=============================================================================


//-----------------------------------------------------------------------------
void Rg2DTri3_::shape(double* H, double r, double s)
{
    H[0] = 1.0 - r - s;
    H[1] = r;
    H[2] = s;
}


//-----------------------------------------------------------------------------
void Rg2DTri3_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[0] = -1; Hs[0] = -1;
    Hr[1] =  1; Hs[1] =  0;
    Hr[2] =  0; Hs[2] =  1;
}

//-----------------------------------------------------------------------------
void Rg2DTri3_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
    Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
    Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
    Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

//=============================================================================
//                          Rg2DTri3G1 
//=============================================================================

//-----------------------------------------------------------------------------
Rg2DTri3G1::Rg2DTri3G1() : Rg2DTri3_(NINT, FE2D_TRI3G1)
{
    const double a = 1.0/3.0;
    gaussPoints[0] = RgGaussPoint(a, a, 0.5);
    init(); 
}

//-----------------------------------------------------------------------------
void Rg2DTri3G1::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[0];
    ao[2] = ai[0];
}

//============================================================================
//                             Rg2DTri6
//============================================================================


//-----------------------------------------------------------------------------
void Rg2DTri6_::shape(double* H, double r, double s)
{
    double r1 = 1.0 - r - s;
    double r2 = r;
    double r3 = s;

    H[0] = r1*(2.0*r1 - 1.0);
    H[1] = r2*(2.0*r2 - 1.0);
    H[2] = r3*(2.0*r3 - 1.0);
    H[3] = 4.0*r1*r2;
    H[4] = 4.0*r2*r3;
    H[5] = 4.0*r3*r1;
}


//-----------------------------------------------------------------------------
void Rg2DTri6_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[0] = -3.0 + 4.0*r + 4.0*s;
    Hr[1] =  4.0*r - 1.0;
    Hr[2] =  0.0;
    Hr[3] =  4.0 - 8.0*r - 4.0*s;
    Hr[4] =  4.0*s;
    Hr[5] = -4.0*s;

    Hs[0] = -3.0 + 4.0*s + 4.0*r;
    Hs[1] =  0.0;
    Hs[2] =  4.0*s - 1.0;
    Hs[3] = -4.0*r;
    Hs[4] =  4.0*r;
    Hs[5] =  4.0 - 8.0*s - 4.0*r;
}

//-----------------------------------------------------------------------------
void Rg2DTri6_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
    Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
    Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
    Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
    Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
    Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
    Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

//=============================================================================
//                          Rg2DTri6G3
//=============================================================================

Rg2DTri6G3::Rg2DTri6G3() : Rg2DTri6_(NINT, FE2D_TRI6G3)
{ 
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 3.0;
    gaussPoints[0] = RgGaussPoint(a, a, a);
    gaussPoints[1] = RgGaussPoint(b, a, a);
    gaussPoints[2] = RgGaussPoint(a, b, a);
    init(); 
}

//-----------------------------------------------------------------------------
void Rg2DTri6G3::project_to_nodes(double* ai, double* ao) const
{
    Matrix H(3, 3);
    for (int n=0; n<3; ++n)
    {
        H(n,0) = 1.0 - gaussPoints[n].getR() - gaussPoints[n].getS();
        H(n,1) = gaussPoints[n].getR();
        H(n,2) = gaussPoints[n].getS();
    }
    H = H.inverse();

    for (int i=0; i<3; ++i)
    {
        ao[i] = 0;
        for (int j=0; j<3; ++j) ao[i] += H(i,j)*ai[j];
    }

    ao[3] = 0.5*(ao[0] + ao[1]);
    ao[4] = 0.5*(ao[1] + ao[2]);
    ao[5] = 0.5*(ao[2] + ao[0]);
}

//=============================================================================
//                          Rg2DQuad4
//=============================================================================


//-----------------------------------------------------------------------------
void Rg2DQuad4_::shape(double* H, double r, double s)
{
    H[0] = 0.25*(1-r)*(1-s);
    H[1] = 0.25*(1+r)*(1-s);
    H[2] = 0.25*(1+r)*(1+s);
    H[3] = 0.25*(1-r)*(1+s);
}


//-----------------------------------------------------------------------------
void Rg2DQuad4_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
    Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
    Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
    Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
}

//-----------------------------------------------------------------------------
void Rg2DQuad4_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
    Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
    Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
    Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
    Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

//=============================================================================
//                          Rg2DQuad4G4 
//=============================================================================

Rg2DQuad4G4::Rg2DQuad4G4() : Rg2DQuad4_(NINT, FE2D_QUAD4G4) 
{
    const double a = 1.0 / sqrt(3.0);
    gaussPoints[0] = RgGaussPoint(-a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, 1);
    init(); 
    m_Ai = m_H.inverse();
}

//-----------------------------------------------------------------------------
void Rg2DQuad4G4::project_to_nodes(double* ai, double* ao) const
{
    int ni = NINT;
    int ne = NELN;
    assert(ni == ne);
    for (int i=0; i<ne; ++i)
    {
        ao[i] = 0;
        for (int j=0; j<ni; ++j) ao[i] += m_Ai[i][j]*ai[j];
    }
}

//=============================================================================
//          Rg2DQuad8 
//=============================================================================


//-----------------------------------------------------------------------------
// shape function at (r,s)
void Rg2DQuad8_::shape(double* H, double r, double s)
{
    H[4] = 0.5*(1 - r*r)*(1 - s);
    H[5] = 0.5*(1 - s*s)*(1 + r);
    H[6] = 0.5*(1 - r*r)*(1 + s);
    H[7] = 0.5*(1 - s*s)*(1 - r);

    H[0] = 0.25*(1 - r)*(1 - s) - 0.5*(H[4] + H[7]);
    H[1] = 0.25*(1 + r)*(1 - s) - 0.5*(H[4] + H[5]);
    H[2] = 0.25*(1 + r)*(1 + s) - 0.5*(H[5] + H[6]);
    H[3] = 0.25*(1 - r)*(1 + s) - 0.5*(H[6] + H[7]);
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
void Rg2DQuad8_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[4] = -r*(1 - s);
    Hr[5] = 0.5*(1 - s*s);
    Hr[6] = -r*(1 + s);
    Hr[7] = -0.5*(1 - s*s);

    Hr[0] = -0.25*(1 - s) - 0.5*(Hr[4] + Hr[7]);
    Hr[1] =  0.25*(1 - s) - 0.5*(Hr[4] + Hr[5]);
    Hr[2] =  0.25*(1 + s) - 0.5*(Hr[5] + Hr[6]);
    Hr[3] = -0.25*(1 + s) - 0.5*(Hr[6] + Hr[7]);

    Hs[4] = -0.5*(1 - r*r);
    Hs[5] = -s*(1 + r);
    Hs[6] = 0.5*(1 - r*r);
    Hs[7] = -s*(1 - r);

    Hs[0] = -0.25*(1 - r) - 0.5*(Hs[4] + Hs[7]);
    Hs[1] = -0.25*(1 + r) - 0.5*(Hs[4] + Hs[5]);
    Hs[2] =  0.25*(1 + r) - 0.5*(Hs[5] + Hs[6]);
    Hs[3] =  0.25*(1 - r) - 0.5*(Hs[6] + Hs[7]);
}


//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
void Rg2DQuad8_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
    Hrr[4] = -(1 - s);
    Hrr[5] = 0.0;
    Hrr[6] = -(1 + s);
    Hrr[7] = 0.0;

    Hrs[4] = r;
    Hrs[5] = -s;
    Hrs[6] = -r;
    Hrs[7] = s;

    Hss[4] = 0.0;
    Hss[5] = -(1 + r);
    Hss[6] = 0.0;
    Hss[7] = -(1 - r);

    Hrr[0] = - 0.5*(Hrr[4] + Hrr[7]);
    Hrr[1] = - 0.5*(Hrr[4] + Hrr[5]);
    Hrr[2] = - 0.5*(Hrr[5] + Hrr[6]);
    Hrr[3] = - 0.5*(Hrr[6] + Hrr[7]);

    Hrs[0] =  0.25 - 0.5*(Hrs[4] + Hrs[7]);
    Hrs[1] = -0.25 - 0.5*(Hrs[4] + Hrs[5]);
    Hrs[2] =  0.25 - 0.5*(Hrs[5] + Hrs[6]);
    Hrs[3] = -0.25 - 0.5*(Hrs[6] + Hrs[7]);

    Hss[0] = - 0.5*(Hss[4] + Hss[7]);
    Hss[1] = - 0.5*(Hss[4] + Hss[5]);
    Hss[2] = - 0.5*(Hss[5] + Hss[6]);
    Hss[3] = - 0.5*(Hss[6] + Hss[7]);
}

//=============================================================================
//       Rg2DQuad8G9
//=============================================================================

Rg2DQuad8G9::Rg2DQuad8G9() : Rg2DQuad8_(NINT, FE2D_QUAD8G9)
{
    // integration point coordinates
    const double a = sqrt(0.6);
    const double w1 = 25.0/81.0;
    const double w2 = 40.0/81.0;
    const double w3 = 64.0/81.0;
    gaussPoints[ 0] = RgGaussPoint(-a, -a, w1);
    gaussPoints[ 1] = RgGaussPoint( 0, -a, w2);
    gaussPoints[ 2] = RgGaussPoint( a, -a, w1);
    gaussPoints[ 3] = RgGaussPoint(-a,  0, w2);
    gaussPoints[ 4] = RgGaussPoint( 0,  0, w3);
    gaussPoints[ 5] = RgGaussPoint( a,  0, w2);
    gaussPoints[ 6] = RgGaussPoint(-a,  a, w1);
    gaussPoints[ 7] = RgGaussPoint( 0,  a, w2);
    gaussPoints[ 8] = RgGaussPoint( a,  a, w1);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN,NELN);
    m_Ai.resize(NELN,NELN);
    A = m_H.transpose()*m_H;
    m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void Rg2DQuad8G9::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i=0; i<NELN; ++i) 
    {
        b[i] = 0;
        for (int j=0; j<NINT; ++j) b[i] += m_H(j,i)*ai[j];
    }

    for (int i=0; i<NELN; ++i) 
    {
        ao[i] = 0;
        for (int j=0; j<NELN; ++j) ao[i] += m_Ai(i,j)*b[j];
    }
}

// Implementation for RgLineElementTraits
RgLineElementTraits::RgLineElementTraits(int ni, int ne, ElementShape es, ElementType et)
    : RgElementTraits(ni, ne, FE_ELEM_TRUSS, es, et)
{
    gaussPoints.resize(ni);

    Gr.resize(ni, ne);
    Grr.resize(ni, ne);
}

void RgLineElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);

    // evaluate shape functions
    std::vector<double> N(m_neln);
    for (int n=0; n<m_nint; ++n)
    {
        shape(N.data(), gaussPoints[n].getR());
        for (int i=0; i<m_neln; ++i) m_H(n, i) = N[i];
    }
    
    // evaluate shape function derivatives
    std::vector<double> Nr(m_neln);
    for (int n=0; n<m_nint; ++n)
    {
        shape_deriv(Nr.data(), gaussPoints[n].getR());
        for (int i=0; i<m_neln; ++i)
        {
            Gr(n, i) = Nr[i];
        }
    }
}

// Implementation for RgLine2_
RgLine2_::RgLine2_(int ni, ElementType et) : RgLineElementTraits(ni, NELN, ET_LINE2, et)
{
}

void RgLine2_::shape(double* H, double r)
{
    H[0] = 0.5*(1.0 - r);
    H[1] = 0.5*(1.0 + r);
}

void RgLine2_::shape_deriv(double* Gr, double r)
{
    Gr[0] = -0.5;
    Gr[1] =  0.5;
}

void RgLine2_::shape_deriv2(double* Grr, double r)
{
    Grr[0] = 0.0;
    Grr[1] = 0.0;
}

// Implementation for RgLine2G1
RgLine2G1::RgLine2G1() : RgLine2_(NINT, FE_LINE2G1)
{
    gaussPoints[0] = RgGaussPoint(0.0, 1.0);
    init();
    
    Matrix H(NELN, NELN);
    for (int i = 0; i < NELN; ++i)
    {
        double r = gaussPoints[0].getR();
        double s = gaussPoints[0].getS();
        auto h = evalH(RgFem::NaturalCoord(r, s));
        for (int j = 0; j < NELN; ++j)
        {
            H(i, j) = h[j];
        }
    }
    
    m_Ai = H.inverse();
}

void RgLine2G1::project_to_nodes(double* ai, double* ao) const
{
    for (int j = 0; j < NELN; ++j)
    {
        ao[j] = 0;
        for (int k = 0; k < NINT; ++k)
        {
            ao[j] += m_Ai[j][k] * ai[k];
        }
    }
}