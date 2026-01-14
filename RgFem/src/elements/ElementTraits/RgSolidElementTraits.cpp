#include "RgSolidElementTraits.h"
#include "RgElementTraits.h"
#include <assert.h>
#include <cmath>
#include "../RgGaussPoint.h"
#include "femcore/RgElementShapeStore.h"
#include "../ElementShape/RgSolidElementShape.h"
#include "elements/NaturalCoord.h"



//=============================================================================
// RgSolidElementTraits
//=============================================================================

RgSolidElementTraits::RgSolidElementTraits(int ni, int ne, ElementShape es, ElementType et) 
    : RgElementTraits(ni, ne, FE_ELEM_SOLID, es, et)
{
    m_shape = nullptr;
    
    // Resize vectors for gauss point data
    gaussPoints.resize(ni);
    
    m_Gr.resize(ni, ne);
    m_Gs.resize(ni, ne);
    m_Gt.resize(ni, ne);
    
    Grr.resize(ni, ne);
    Gsr.resize(ni, ne);
    Gtr.resize(ni, ne);
        
    Grs.resize(ni, ne);
    Gss.resize(ni, ne);
    Gts.resize(ni, ne);
        
    Grt.resize(ni, ne);
    Gst.resize(ni, ne);
    Gtt.resize(ni, ne);
    
    // Set number of faces based on element shape
    m_faces = 0;
    switch (es)
    {
    case ET_TET4:
    case ET_TET10:
        m_faces = 4;
        break;
    case ET_PENTA6:
    case ET_PENTA15:
        m_faces = 5;
        break;
    case ET_HEX8:
    case ET_HEX20:
    case ET_HEX27:
        m_faces = 6;
        break;
    case ET_PYRA5:
    case ET_PYRA13:
        m_faces = 5;
        break;
    default:
        assert(false);
    }
}

void RgSolidElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);
    const int MAX_NODES = 27; // Maximum nodes for HEX27
    
    // get shape class
    m_shape = dynamic_cast<RgSolidElementShape*>(RgElementShapeStore::GetInstance()->GetElementShape(m_spec.eshape));
    
    // calculate shape function values at gauss points
    std::vector<double> N(MAX_NODES,0);
    for (int n=0; n<m_nint; ++n)
    {
        NaturalCoord coord(gaussPoints[n]);
        N = m_shape->evalH(coord);
        for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
    }
    
    // calculate local derivatives of shape functions at gauss points
    for (int n=0; n<m_nint; ++n)
    {
        NaturalCoord coord(gaussPoints[n]);
        auto result = m_shape->evalDeriv(coord); // result is {dN/dr, dN/ds, dN/dt}
        for (int i=0; i<m_neln; ++i)
        {
            m_Gr[n][i] = result[0][i];
            m_Gs[n][i] = result[1][i];
            m_Gt[n][i] = result[2][i];
        }
    }
    
    // calculate local second derivatives of shape functions at gauss points
    double Hrr[MAX_NODES], Hss[MAX_NODES], Htt[MAX_NODES], Hrs[MAX_NODES], Hst[MAX_NODES], Hrt[MAX_NODES];
    for (int n=0; n<m_nint; ++n)
    {
        NaturalCoord coord(gaussPoints[n]);
        auto result = m_shape->evalDeriv2(coord); // 6  Nκ,  ע˳
        for (int i=0; i<m_neln; ++i)
        {
            Grr[n][i] = result[0][i]; 
            Gss[n][i] = result[1][i];
            Gtt[n][i] = result[2][i];

            Grs[n][i] = result[3][i];
            Gsr[n][i] = result[3][i];

            Grt[n][i] = result[4][i];
            Gtr[n][i] = result[4][i];
            
            Gst[n][i] = result[5][i];
            Gts[n][i] = result[5][i];
            
        }
    }
}

std::vector<double> RgSolidElementTraits::evalH(const NaturalCoord& coord)
{
    return m_shape->evalH(coord);
}

std::vector<std::vector<double>> RgSolidElementTraits::evalDeriv(const NaturalCoord& coord)
{
    return m_shape->evalDeriv(coord);
}

std::vector<std::vector<double>> RgSolidElementTraits::evalDeriv2(const NaturalCoord& coord)
{
    return m_shape->evalDeriv2(coord);
}

void RgSolidElementTraits::project_to_nodes(double* ai, double* ao) const
{
    // Default implementation is empty
}

const RgGaussPoint RgSolidElementTraits::gaussPoint(int n)
{
    if (n >= 0 && n < static_cast<int>(gaussPoints.size())) {
        return gaussPoints[n];
    }
    // Return a default constructed RgGaussPoint if index is out of bounds
    return RgGaussPoint();
}

//void RgSolidElementTraits::shape_fnc(double* H, double r, double s, double t)
//{
//    m_shape->shape_fnc(H, r, s, t); 
//}
//
//void RgSolidElementTraits::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
//{
//    m_shape->shape_deriv(Hr, Hs, Ht, r, s, t);
//}
//
//void RgSolidElementTraits::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
//{
//    m_shape->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t);
//}

//=============================================================================
//! Base class for 8-node hexahedral elements
//=============================================================================

RgHex8_::RgHex8_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_HEX8, et) 
{
}

void RgHex8_::init()
{
    // TODO: Implementation for shape classes allocation if needed
}

//=============================================================================
// 8-node hexahedral elements with 8-point gaussian quadrature
//=============================================================================

RgHex8G8::RgHex8G8() : RgHex8_(NINT, FE_HEX8G8) 
{
    // integration point coordinates
    const double a = 1.0 / sqrt(3.0);
    gaussPoints[0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[4] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[5] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[6] = RgGaussPoint( a,  a,  a, 1);
    gaussPoints[7] = RgGaussPoint(-a,  a,  a, 1);
    
    init();
    m_Hi = m_H.inverse();
}

void RgHex8G8::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k) 
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}

//=============================================================================
// 8-node hexahedral elements with 6-point reduced integration rule
//=============================================================================

RgHex8RI::RgHex8RI() : RgHex8_(NINT, FE_HEX8RI) 
{
    // This is for a six point integration rule
    // integration point coordinates
    const double a = 8.0 / 6.0;
    gaussPoints[0] = RgGaussPoint(-1,  0,  0, a);
    gaussPoints[1] = RgGaussPoint( 1,  0,  0, a);
    gaussPoints[2] = RgGaussPoint( 0, -1,  0, a);
    gaussPoints[3] = RgGaussPoint( 0,  1,  0, a);
    gaussPoints[4] = RgGaussPoint( 0,  0, -1, a);
    gaussPoints[5] = RgGaussPoint( 0,  0,  1, a);
    
    init();
}

void RgHex8RI::project_to_nodes(double* ai, double* ao) const
{
    // TODO: implement this
}

//=============================================================================
// 8-node hexahedral element with uniform deformation gradient
//=============================================================================

RgHex8G1::RgHex8G1() : RgHex8_(NINT, FE_HEX8G1) 
{
    // single gauss-point integration rule
    gaussPoints[0] = RgGaussPoint(0, 0, 0, 8.0);
    
    init();
}

void RgHex8G1::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[0];
    ao[2] = ai[0];
    ao[3] = ai[0];
    ao[4] = ai[0];
    ao[5] = ai[0];
    ao[6] = ai[0];
    ao[7] = ai[0];
}

//=============================================================================
//! Base class for 4-node linear tetrahedrons
//=============================================================================

RgTet4_::RgTet4_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_TET4, et) 
{
}

void RgTet4_::init()
{
    // TODO: Implementation for shape classes allocation if needed
}

//=============================================================================
// single Gauss point integrated tet element
//=============================================================================

RgTet4G1::RgTet4G1() : RgTet4_(NINT, FE_TET4G1) 
{
    // gaussian integration for tetrahedral elements
    const double a = 0.25;
    const double w = 1.0 / 6.0;
    
    gaussPoints[0] = RgGaussPoint(a, a, a, w);
    
    init();
}

void RgTet4G1::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[0];
    ao[2] = ai[0];
    ao[3] = ai[0];
}

//=============================================================================
// 4-node tetrahedral element using a 4-node Gaussian integration rule
//=============================================================================

RgTet4G4::RgTet4G4() : RgTet4_(NINT, FE_TET4G4) 
{
    // gaussian integration for tetrahedral elements
    const double a = 0.58541020;
    const double b = 0.13819660;
    const double w = 1.0 / 24.0;
    
    gaussPoints[0] = RgGaussPoint(b, b, b, w);
    gaussPoints[1] = RgGaussPoint(a, b, b, w);
    gaussPoints[2] = RgGaussPoint(b, a, b, w);
    gaussPoints[3] = RgGaussPoint(b, b, a, w);
    
    init();
    m_Hi = m_H.inverse();
}

void RgTet4G4::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k) 
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}


//=============================================================================
//                        P E N T A 6 G 6
//=============================================================================

RgPenta6_::RgPenta6_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_PENTA6, et)
{
}

void RgPenta6_::init()
{
    // TODO: Implementation for shape classes allocation if needed
}

RgPenta6G6::RgPenta6G6() : RgPenta6_(NINT, FE_PENTA6G6)
{
    // gaussian integration for pentahedral elements
    const double a = 0.577350269189626;
    const double w = 1.0 / 6.0;
    
    gaussPoints[0] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[1] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[2] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[3] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[4] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[5] = RgGaussPoint(0, 0.5,  a, w);

    init();
    m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void RgPenta6G6::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k) 
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}

//=============================================================================
//                        P E N T A 1 5 G 2 1
//=============================================================================

RgPenta15_::RgPenta15_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_PENTA15, et)
{
}

void RgPenta15_::init()
{
    // TODO: Implementation for shape classes allocation if needed
}

RgPenta15G21::RgPenta15G21() : RgPenta15_(NINT, FE_PENTA15G21)
{
    // gaussian integration for pentahedral elements
    const double a = 0.577350269189626;
    const double w = 1.0 / 6.0;
    
    gaussPoints[ 0] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[ 1] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[ 2] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[ 3] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[ 4] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[ 5] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[ 6] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[ 7] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[ 8] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[ 9] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[10] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[11] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[12] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[13] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[14] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[15] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[16] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[17] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[18] = RgGaussPoint(0, 0.5, -a, w);
    gaussPoints[19] = RgGaussPoint(0, 0.5,  a, w);
    gaussPoints[20] = RgGaussPoint(0, 0.5, -a, w);

    init();
    m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void RgPenta15G21::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k) 
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}

//=============================================================================
//              H E X 2 0
//=============================================================================

int RgHex20_::Nodes(int order)
{
    switch (order)
    {
    case 2: return 20; break;
    case 1: return 8; break;
    case 0: return 1; break;
    default:
        assert(false);
        return 20;
    }
}

void RgHex20_::init()
{
    // TODO: Implementation for shape classes allocation if needed
}

//=============================================================================
//              H E X 2 0 G 8
//=============================================================================

int RgHex20G8::ni[NELN] = {};

RgHex20G8::RgHex20G8() : RgHex20_(NINT, FE_HEX20G8)
{
    // integration point coordinates
    const double a = 1.0 / sqrt(3.0);
    gaussPoints[0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[4] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[5] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[6] = RgGaussPoint( a,  a,  a, 1);
    gaussPoints[7] = RgGaussPoint(-a,  a,  a, 1);
    
    init();
    
    MT.resize(NELN, NINT);
    for (int i=0; i<NINT; ++i)
        for (int n=0; n<NELN; ++n)
            MT(n,i) = m_H(i,n);
    
    Hi.resize(NELN, NELN);
    Hi = MT*MT.transpose();
    Hi = Hi.inverse();
}

//-----------------------------------------------------------------------------
//! Use least-squares extrapolation
void RgHex20G8::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i=0; i<NELN; ++i)
    {
        b[i] = 0;
        for (int j=0; j<NINT; ++j) b[i] += MT(i,j)*ai[j];
    }

    for (int i=0; i<NELN; ++i)
    {
        ao[i] = 0.0;
        for (int j=0; j<NELN; ++j) ao[i] += Hi[i][j]*b[j];
    }
}

//=============================================================================
//              H E X 2 0 G 2 7
//=============================================================================

int RgHex20G27::ni[NELN] = {};

RgHex20G27::RgHex20G27() : RgHex20_(NINT, FE_HEX20G27)
{
    // integration point coordinates
    const double a = sqrt(0.6);
    gaussPoints[ 0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[ 1] = RgGaussPoint( 0, -a, -a, 1);
    gaussPoints[ 2] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[ 3] = RgGaussPoint(-a,  0, -a, 1);
    gaussPoints[ 4] = RgGaussPoint( 0,  0, -a, 1);
    gaussPoints[ 5] = RgGaussPoint( a,  0, -a, 1);
    gaussPoints[ 6] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[ 7] = RgGaussPoint( 0,  a, -a, 1);
    gaussPoints[ 8] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[ 9] = RgGaussPoint(-a, -a,  0, 1);
    gaussPoints[10] = RgGaussPoint( 0, -a,  0, 1);
    gaussPoints[11] = RgGaussPoint( a, -a,  0, 1);
    gaussPoints[12] = RgGaussPoint(-a,  0,  0, 1);
    gaussPoints[13] = RgGaussPoint( 0,  0,  0, 1);
    gaussPoints[14] = RgGaussPoint( a,  0,  0, 1);
    gaussPoints[15] = RgGaussPoint(-a,  a,  0, 1);
    gaussPoints[16] = RgGaussPoint( 0,  a,  0, 1);
    gaussPoints[17] = RgGaussPoint( a,  a,  0, 1);
    gaussPoints[18] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[19] = RgGaussPoint( 0, -a,  a, 1);
    gaussPoints[20] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[21] = RgGaussPoint(-a,  0,  a, 1);
    gaussPoints[22] = RgGaussPoint( 0,  0,  a, 1);
    gaussPoints[23] = RgGaussPoint( a,  0,  a, 1);
    gaussPoints[24] = RgGaussPoint(-a,  a,  a, 1);
    gaussPoints[25] = RgGaussPoint( 0,  a,  a, 1);
    gaussPoints[26] = RgGaussPoint( a,  a,  a, 1);

    init();

    MT.resize(NELN, NINT);
    for (int i=0; i<NINT; ++i)
        for (int n=0; n<NELN; ++n)
            MT(n,i) = m_H(i,n);

    Hi.resize(NELN, NELN);
    Hi = MT*MT.transpose();
    Hi = Hi.inverse();
}

//-----------------------------------------------------------------------------
//! Use least-squares extrapolation
void RgHex20G27::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i=0; i<NELN; ++i)
    {
        b[i] = 0;
        for (int j=0; j<NINT; ++j) b[i] += MT(i,j)*ai[j];
    }

    for (int i=0; i<NELN; ++i)
    {
        ao[i] = 0.0;
        for (int j=0; j<NELN; ++j) ao[i] += Hi[i][j]*b[j];
    }
}

//=============================================================================
//              H E X 2 7 G 2 7
//=============================================================================

RgHex27_::RgHex27_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_HEX27, et)
{
}

RgHex27G27::RgHex27G27() : RgHex27_(NINT, FE_HEX27G27)
{
    // integration point coordinates
    const double a = sqrt(0.6);
    gaussPoints[ 0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[ 1] = RgGaussPoint( 0, -a, -a, 1);
    gaussPoints[ 2] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[ 3] = RgGaussPoint(-a,  0, -a, 1);
    gaussPoints[ 4] = RgGaussPoint( 0,  0, -a, 1);
    gaussPoints[ 5] = RgGaussPoint( a,  0, -a, 1);
    gaussPoints[ 6] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[ 7] = RgGaussPoint( 0,  a, -a, 1);
    gaussPoints[ 8] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[ 9] = RgGaussPoint(-a, -a,  0, 1);
    gaussPoints[10] = RgGaussPoint( 0, -a,  0, 1);
    gaussPoints[11] = RgGaussPoint( a, -a,  0, 1);
    gaussPoints[12] = RgGaussPoint(-a,  0,  0, 1);
    gaussPoints[13] = RgGaussPoint( 0,  0,  0, 1);
    gaussPoints[14] = RgGaussPoint( a,  0,  0, 1);
    gaussPoints[15] = RgGaussPoint(-a,  a,  0, 1);
    gaussPoints[16] = RgGaussPoint( 0,  a,  0, 1);
    gaussPoints[17] = RgGaussPoint( a,  a,  0, 1);
    gaussPoints[18] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[19] = RgGaussPoint( 0, -a,  a, 1);
    gaussPoints[20] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[21] = RgGaussPoint(-a,  0,  a, 1);
    gaussPoints[22] = RgGaussPoint( 0,  0,  a, 1);
    gaussPoints[23] = RgGaussPoint( a,  0,  a, 1);
    gaussPoints[24] = RgGaussPoint(-a,  a,  a, 1);
    gaussPoints[25] = RgGaussPoint( 0,  a,  a, 1);
    gaussPoints[26] = RgGaussPoint( a,  a,  a, 1);

    init();
    m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void RgHex27G27::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k) 
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}

//=============================================================================
//              P Y R A 5 G 8
//=============================================================================

RgPyra5_::RgPyra5_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_PYRA5, et)
{
}

RgPyra5G8::RgPyra5G8() : RgPyra5_(NINT, FE_PYRA5G8)
{
    // integration point coordinates
    const double a = 0.577350269189626;
    gaussPoints[0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[4] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[5] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[6] = RgGaussPoint( a,  a,  a, 1);
    gaussPoints[7] = RgGaussPoint(-a,  a,  a, 1);

    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN,NELN);
    m_Ai.resize(NELN,NELN);
    A = m_H.transpose()*m_H;
    m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void RgPyra5G8::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(5);
    for (int i=0; i<5; ++i)
    {
        b[i] = 0;
        for (int j=0; j<NINT; ++j) b[i] += m_H(j,i)*ai[j];
    }

    for (int i=0; i<5; ++i)
    {
        ao[i] = 0.0;
        for (int j=0; j<5; ++j) ao[i] += m_Ai[i][j]*b[j];
    }
}

//=============================================================================
//              P Y R A 1 3 G 8
//=============================================================================

RgPyra13_::RgPyra13_(int ni, ElementType et) : RgSolidElementTraits(ni, NELN, ET_PYRA13, et)
{
}

RgPyra13G8::RgPyra13G8() : RgPyra13_(NINT, FE_PYRA13G8)
{
    // integration point coordinates
    const double a = 0.577350269189626;
    gaussPoints[0] = RgGaussPoint(-a, -a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, -a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, -a, 1);
    gaussPoints[4] = RgGaussPoint(-a, -a,  a, 1);
    gaussPoints[5] = RgGaussPoint( a, -a,  a, 1);
    gaussPoints[6] = RgGaussPoint( a,  a,  a, 1);
    gaussPoints[7] = RgGaussPoint(-a,  a,  a, 1);
    
    init();
    
    // we need Ai to project integration point data to the nodes
    Matrix A(NELN,NELN);
    m_Ai.resize(NELN,NELN);
    A = m_H.transpose()*m_H;
    m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void RgPyra13G8::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(13);
    for (int i=0; i<13; ++i)
    {
        b[i] = 0;
        for (int j=0; j<NINT; ++j) b[i] += m_H(j,i)*ai[j];
    }
    
    for (int i=0; i<13; ++i)
    {
        ao[i] = 0.0;
        for (int j=0; j<13; ++j) ao[i] += m_Ai[i][j]*b[j];
    }
}