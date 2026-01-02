#include "RgSurfaceElementTraits.h"
#include "RgElementTraits.h"
#include <assert.h>
#include <cmath>
#include "../RgGaussPoint.h"
#include "femcore/RgElementShapeStore.h"
#include "elements/ElementShape/RgSurfaceElementShape.h"
#include "elements/NaturalCoord.h"

using namespace RgFem;

// RgSurfaceElementTraits implementation
RgSurfaceElementTraits::RgSurfaceElementTraits(int ni, int ne, ElementShape es, ElementType et)
    : RgElementTraits(ni, ne, FE_ELEM_SURFACE, es, et)
{
    m_shape = nullptr;
    
    // Resize vectors for gauss point data
    gaussPoints.resize(ni);
    
    m_Gr.resize(ni, ne);
    m_Gs.resize(ni, ne);
    
    Grr.resize(ni, ne);
    Gsr.resize(ni, ne);
    Grs.resize(ni, ne);
    Gss.resize(ni, ne);
    
    // Initialize parametric coordinates of element center
    cr = 0.0;
    cs = 0.0;
}

std::vector<double> RgSurfaceElementTraits::evalH(const RgFem::NaturalCoord& coord)
{
    return m_shape->evalH(coord);
}

std::vector<std::vector<double>> RgSurfaceElementTraits::evalDeriv(const RgFem::NaturalCoord& coord)
{
    return m_shape->evalDeriv(coord);
}

std::vector<std::vector<double>> RgSurfaceElementTraits::evalDeriv2(const RgFem::NaturalCoord& coord)
{
    return m_shape->evalDeriv2(coord);
}

void RgSurfaceElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);
    const int MAX_NODES = 27; // Maximum nodes for high-order elements
    
    // get shape class
    m_shape = dynamic_cast<RgSurfaceElementShape*>(RgElementShapeStore::GetInstance()->GetElementShape(m_spec.eshape));
    assert(m_shape && (m_shape->shapeType() == m_spec.eshape));
    
    // calculate shape function values at gauss points
    std::vector<double> N(MAX_NODES, 0);
    for (int n=0; n<m_nint; ++n)
    {
        RgFem::NaturalCoord coord(gaussPoints[n]);
        N = m_shape->evalH(coord);
        for (int i=0; i<m_neln; ++i) m_H(n, i) = N[i];
    }
    
    // calculate local derivatives of shape functions at gauss points
    for (int n=0; n<m_nint; ++n)
    {
        RgFem::NaturalCoord coord(gaussPoints[n]);
        auto result = m_shape->evalDeriv(coord); // result is {dN/dr, dN/ds}
        for (int i=0; i<m_neln; ++i)
        {
            m_Gr(n, i) = result[0][i];  // dr derivatives
            m_Gs(n, i) = result[1][i];  // ds derivatives
        }
    }
    
    // calculate local second derivatives of shape functions at gauss points
    for (int n=0; n<m_nint; ++n)
    {
        RgFem::NaturalCoord coord(gaussPoints[n]);
        auto result = m_shape->evalDeriv2(coord); // result is {d2N/dr2, d2N/ds2, d2N/drdt}
        for (int i=0; i<m_neln; ++i)
        {
            Grr(n, i) = result[0][i];  // d2N/dr2
            Gss(n, i) = result[1][i];  // d2N/ds2
            Grs(n, i) = result[2][i];  // d2N/drdt
            Gsr(n, i) = result[2][i];  // d2N/dsdr (symmetric to Grs)
        }
    }
}

// Implementation for RgQuad4_
void RgQuad4_::init()
{
    // allocate shape classes if needed
    // centroid coordinates
    cr = 0.0;
    cs = 0.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgQuad4G4
RgQuad4G4::RgQuad4G4() : RgQuad4_(NINT, FE_QUAD4G4)
{
    const double a = 1.0 / sqrt(3.0);
    gaussPoints[0] = RgGaussPoint(-a, -a, 1);
    gaussPoints[1] = RgGaussPoint( a, -a, 1);
    gaussPoints[2] = RgGaussPoint( a,  a, 1);
    gaussPoints[3] = RgGaussPoint(-a,  a, 1);
    init();
    m_Hi = m_H.inverse();
}

void RgQuad4G4::project_to_nodes(double* ai, double* ao) const
{
    int ni = NINT;
    int ne = NELN;
    assert(ni == ne);
    for (int i=0; i<ne; ++i)
    {
        ao[i] = 0;
        for (int j=0; j<ni; ++j) ao[i] += m_Hi[i][j]*ai[j];
    }
}

// Implementation for RgQuad4NI
RgQuad4NI::RgQuad4NI() : RgQuad4_(NINT, FE_QUAD4NI)
{
    gaussPoints[0] = RgGaussPoint(-1, -1, 1);
    gaussPoints[1] = RgGaussPoint( 1, -1, 1);
    gaussPoints[2] = RgGaussPoint( 1,  1, 1);
    gaussPoints[3] = RgGaussPoint(-1,  1, 1);
    init();
}

void RgQuad4NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
}

// Implementation for RgTri3_
void RgTri3_::init()
{
    // centroid coordinates
    cr = 1.0 / 3.0;
    cs = 1.0 / 3.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgTri3G1
RgTri3G1::RgTri3G1() : RgTri3_(NINT, FE_TRI3G1)
{
    const double a = 1.0 / 3.0;
    gaussPoints[0] = RgGaussPoint(a, a, 0.5);
    init();
}

void RgTri3G1::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[0];
    ao[2] = ai[0];
}

// Implementation for RgTri3G3
RgTri3G3::RgTri3G3() : RgTri3_(NINT, FE_TRI3G3)
{
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 3.0;
    gaussPoints[0] = RgGaussPoint(a, a, a);
    gaussPoints[1] = RgGaussPoint(b, a, a);
    gaussPoints[2] = RgGaussPoint(a, b, a);
    init();
    m_Hi = m_H.inverse();
}

void RgTri3G3::project_to_nodes(double* ai, double* ao) const
{
    int ni = NINT;
    int ne = NELN;
    assert(ni == ne);
    for (int i=0; i<ne; ++i)
    {
        ao[i] = 0;
        for (int j=0; j<ni; ++j) ao[i] += m_Hi[i][j]*ai[j];
    }
}

// Implementation for RgTri3G7
RgTri3G7::RgTri3G7() : RgTri3_(NINT, FE_TRI3G7)
{
    const double w = 1.0 / 2.0;
    gaussPoints[0] = RgGaussPoint(0.333333333333333, 0.333333333333333, w * 0.225000000000000);
    gaussPoints[1] = RgGaussPoint(0.797426985353087, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[2] = RgGaussPoint(0.101286507323456, 0.797426985353087, w * 0.125939180544827);
    gaussPoints[3] = RgGaussPoint(0.101286507323456, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[4] = RgGaussPoint(0.470142064105115, 0.470142064105115, w * 0.132394152788506);
    gaussPoints[5] = RgGaussPoint(0.470142064105115, 0.059715871789770, w * 0.132394152788506);
    gaussPoints[6] = RgGaussPoint(0.059715871789770, 0.470142064105115, w * 0.132394152788506);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgTri3G7::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i = 0; i < NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j < NINT; ++j) b[i] += m_H(j, i) * ai[j];
    }

    for (int i = 0; i < NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < NELN; ++j) ao[i] += m_Ai(i, j) * b[j];
    }
}

// Implementation for RgTri3NI
RgTri3NI::RgTri3NI() : RgTri3_(NINT, FE_TRI3NI)
{
    gaussPoints[0] = RgGaussPoint(0, 0, 0.0);
    gaussPoints[1] = RgGaussPoint(1, 0, 0.0);
    gaussPoints[2] = RgGaussPoint(0, 1, 0.0);
    init();
}

void RgTri3NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
}

// Implementation for RgTri6_
void RgTri6_::init()
{
    // centroid coordinates
    cr = 1.0 / 3.0;
    cs = 1.0 / 3.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgTri6G3
RgTri6G3::RgTri6G3() : RgTri6_(NINT, FE_TRI6G3)
{
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 3.0;
    gaussPoints[0] = RgGaussPoint(a, a, a);
    gaussPoints[1] = RgGaussPoint(b, a, a);
    gaussPoints[2] = RgGaussPoint(a, b, a);
    init();
}

void RgTri6G3::project_to_nodes(double* ai, double* ao) const
{
    Matrix H(3, 3);
    for (int n = 0; n < 3; ++n)
    {
        H(n, 0) = 1.0 - gaussPoints[n].getR() - gaussPoints[n].getS();
        H(n, 1) = gaussPoints[n].getR();
        H(n, 2) = gaussPoints[n].getS();
    }
    H = H.inverse();

    for (int i = 0; i < 3; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < 3; ++j) ao[i] += H(i, j) * ai[j];
    }

    ao[3] = 0.5 * (ao[0] + ao[1]);
    ao[4] = 0.5 * (ao[1] + ao[2]);
    ao[5] = 0.5 * (ao[2] + ao[0]);
}

// Implementation for RgTri6G4
RgTri6G4::RgTri6G4() : RgTri6_(NINT, FE_TRI6G4)
{
    const double a = 1.0 / 3.0;
    const double b = 1.0 / 5.0;
    const double c = 3.0 / 5.0;
    gaussPoints[0] = RgGaussPoint(a, a, -27.0 / 96.0);
    gaussPoints[1] = RgGaussPoint(c, b, 25.0 / 96.0);
    gaussPoints[2] = RgGaussPoint(b, b, 25.0 / 96.0);
    gaussPoints[3] = RgGaussPoint(b, c, 25.0 / 96.0);
    init();
}

void RgTri6G4::project_to_nodes(double* ai, double* ao) const
{
    // TODO: Implement this
}

// Implementation for RgTri6G7
RgTri6G7::RgTri6G7() : RgTri6_(NINT, FE_TRI6G7)
{
    const double w = 1.0 / 2.0;
    gaussPoints[0] = RgGaussPoint(0.333333333333333, 0.333333333333333, w * 0.225000000000000);
    gaussPoints[1] = RgGaussPoint(0.797426985353087, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[2] = RgGaussPoint(0.101286507323456, 0.797426985353087, w * 0.125939180544827);
    gaussPoints[3] = RgGaussPoint(0.101286507323456, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[4] = RgGaussPoint(0.470142064105115, 0.470142064105115, w * 0.132394152788506);
    gaussPoints[5] = RgGaussPoint(0.470142064105115, 0.059715871789770, w * 0.132394152788506);
    gaussPoints[6] = RgGaussPoint(0.059715871789770, 0.470142064105115, w * 0.132394152788506);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgTri6G7::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i = 0; i < NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j < NINT; ++j) b[i] += m_H(j, i) * ai[j];
    }

    for (int i = 0; i < NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < NELN; ++j) ao[i] += m_Ai(i, j) * b[j];
    }
}

// Implementation for RgTri6NI
RgTri6NI::RgTri6NI() : RgTri6_(NINT, FE_TRI6NI)
{
	const double a = 0.0;
	const double b = 1.0 / 6.0;
	gaussPoints[0] = RgGaussPoint(0.0, 0.0, a);
	gaussPoints[1] = RgGaussPoint(1.0, 0.0, a);
	gaussPoints[2] = RgGaussPoint(0.0, 1.0, a);
	gaussPoints[3] = RgGaussPoint(0.5, 0.0, b);
	gaussPoints[4] = RgGaussPoint(0.5, 0.5, b);
	gaussPoints[5] = RgGaussPoint(0.0, 0.5, b);
	init();
}

void RgTri6NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
    ao[4] = ai[4];
    ao[5] = ai[5];
}

// Implementation for RgTri7_
void RgTri7_::init()
{
    // centroid coordinates
    cr = 1.0 / 3.0;
    cs = 1.0 / 3.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgTri7G3
RgTri7G3::RgTri7G3() : RgTri7_(NINT, FE_TRI7G3)
{
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 3.0;
    gaussPoints[0] = RgGaussPoint(a, a, a);
    gaussPoints[1] = RgGaussPoint(b, a, a);
    gaussPoints[2] = RgGaussPoint(a, b, a);
    init();
}

void RgTri7G3::project_to_nodes(double* ai, double* ao) const
{
	Matrix H(3, 3);
	for (int n = 0; n < 3; ++n)
	{
		H(n, 0) = 1.0 - gaussPoints[n].getR() - gaussPoints[n].getS();
		H(n, 1) = gaussPoints[n].getR();
		H(n, 2) = gaussPoints[n].getS();
	}
	H.inverse();

	for (int i = 0; i < 3; ++i)
	{
		ao[i] = 0;
		for (int j = 0; j < 3; ++j) ao[i] += H(i, j) * ai[j];
	}

	ao[3] = 0.5 * (ao[0] + ao[1]);
	ao[4] = 0.5 * (ao[1] + ao[2]);
	ao[5] = 0.5 * (ao[2] + ao[0]);
	ao[6] = (ao[0] + ao[1] + ao[2]) / 3.0;
}

// Implementation for RgTri7G4
RgTri7G4::RgTri7G4() : RgTri7_(NINT, FE_TRI7G4)
{
    const double a = 1.0 / 3.0;
    const double b = 1.0 / 5.0;
    const double c = 3.0 / 5.0;
    gaussPoints[0] = RgGaussPoint(a, a, -27.0 / 96.0);
    gaussPoints[1] = RgGaussPoint(c, b, 25.0 / 96.0);
    gaussPoints[2] = RgGaussPoint(b, b, 25.0 / 96.0);
    gaussPoints[3] = RgGaussPoint(b, c, 25.0 / 96.0);
    init();
}

void RgTri7G4::project_to_nodes(double* ai, double* ao) const
{
    // TODO: Implement this
}

// Implementation for RgTri7G7
RgTri7G7::RgTri7G7() : RgTri7_(NINT, FE_TRI7G7)
{
    const double w = 1.0 / 2.0;
    gaussPoints[0] = RgGaussPoint(0.333333333333333, 0.333333333333333, w * 0.225000000000000);
    gaussPoints[1] = RgGaussPoint(0.797426985353087, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[2] = RgGaussPoint(0.101286507323456, 0.797426985353087, w * 0.125939180544827);
    gaussPoints[3] = RgGaussPoint(0.101286507323456, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[4] = RgGaussPoint(0.470142064105115, 0.470142064105115, w * 0.132394152788506);
    gaussPoints[5] = RgGaussPoint(0.470142064105115, 0.059715871789770, w * 0.132394152788506);
    gaussPoints[6] = RgGaussPoint(0.059715871789770, 0.470142064105115, w * 0.132394152788506);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgTri7G7::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i = 0; i < NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j < NINT; ++j) b[i] += m_H(j, i) * ai[j];
    }

    for (int i = 0; i < NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < NELN; ++j) ao[i] += m_Ai(i, j) * b[j];
    }
}

// Implementation for RgTri10_
void RgTri10_::init()
{
    // centroid coordinates
    cr = 1.0 / 3.0;
    cs = 1.0 / 3.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgTri10G7
RgTri10G7::RgTri10G7() : RgTri10_(NINT, FE_TRI10G7)
{
    const double w = 1.0 / 2.0;
    gaussPoints[0] = RgGaussPoint(0.333333333333333, 0.333333333333333, w * 0.225000000000000);
    gaussPoints[1] = RgGaussPoint(0.797426985353087, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[2] = RgGaussPoint(0.101286507323456, 0.797426985353087, w * 0.125939180544827);
    gaussPoints[3] = RgGaussPoint(0.101286507323456, 0.101286507323456, w * 0.125939180544827);
    gaussPoints[4] = RgGaussPoint(0.470142064105115, 0.470142064105115, w * 0.132394152788506);
    gaussPoints[5] = RgGaussPoint(0.470142064105115, 0.059715871789770, w * 0.132394152788506);
    gaussPoints[6] = RgGaussPoint(0.059715871789770, 0.470142064105115, w * 0.132394152788506);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgTri10G7::project_to_nodes(double* ai, double* ao) const
{
    // TODO: Implement this
}

// Implementation for RgTri10G12
RgTri10G12::RgTri10G12() : RgTri10_(NINT, FE_TRI10G12)
{
    gaussPoints[0] = RgGaussPoint(0.063089014, 0.873821971, 0.025422453);
    gaussPoints[1] = RgGaussPoint(0.873821971, 0.063089014, 0.025422453);
    gaussPoints[2] = RgGaussPoint(0.063089014, 0.063089014, 0.025422453);
    gaussPoints[3] = RgGaussPoint(0.249286745, 0.501426510, 0.058393138);
    gaussPoints[4] = RgGaussPoint(0.501426510, 0.249286745, 0.058393138);
    gaussPoints[5] = RgGaussPoint(0.249286745, 0.249286745, 0.058393138);
    gaussPoints[6] = RgGaussPoint(0.053145050, 0.636502499, 0.041425538);
    gaussPoints[7] = RgGaussPoint(0.636502499, 0.053145050, 0.041425538);
    gaussPoints[8] = RgGaussPoint(0.310352451, 0.636502499, 0.041425538);
    gaussPoints[9] = RgGaussPoint(0.636502499, 0.310352451, 0.041425538);
    gaussPoints[10] = RgGaussPoint(0.310352451, 0.053145050, 0.041425538);
    gaussPoints[11] = RgGaussPoint(0.053145050, 0.310352451, 0.041425538);
    init();

    // we need Ai to project integration point data to the nodes
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgTri10G12::project_to_nodes(double* ai, double* ao) const
{
    // TODO: Implement this
}

// Implementation for RgQuad8_
void RgQuad8_::init()
{
    // centroid coordinates
    cr = 0.0;
    cs = 0.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgQuad8G9
RgQuad8G9::RgQuad8G9() : RgQuad8_(NINT, FE_QUAD8G9)
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
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgQuad8G9::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i = 0; i < NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j < NINT; ++j) b[i] += m_H(j, i) * ai[j];
    }

    for (int i = 0; i < NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < NELN; ++j) ao[i] += m_Ai(i, j) * b[j];
    }
}

// Implementation for RgQuad8NI
RgQuad8NI::RgQuad8NI() : RgQuad8_(NINT, FE_QUAD8NI)
{
    double w = 1. / 9.;
    gaussPoints[0] = RgGaussPoint(-1, -1, w);
    gaussPoints[1] = RgGaussPoint( 1, -1, w);
    gaussPoints[2] = RgGaussPoint( 1,  1, w);
    gaussPoints[3] = RgGaussPoint(-1,  1, w);
    gaussPoints[4] = RgGaussPoint( 0, -1, 4 * w);
    gaussPoints[5] = RgGaussPoint( 1,  0, 4 * w);
    gaussPoints[6] = RgGaussPoint( 0,  1, 4 * w);
    gaussPoints[7] = RgGaussPoint(-1,  0, 4 * w);
    init();
}

void RgQuad8NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
    ao[4] = ai[4];
    ao[5] = ai[5];
    ao[6] = ai[6];
    ao[7] = ai[7];
}

// Implementation for RgQuad9_
void RgQuad9_::init()
{
    // centroid coordinates
    cr = 0.0;
    cs = 0.0;

    // initialize base class
    RgSurfaceElementTraits::init();
}

// Implementation for RgQuad9G9
RgQuad9G9::RgQuad9G9() : RgQuad9_(NINT, FE_QUAD9G9)
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
    Matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose() * m_H;
    m_Ai = A.inverse();
}

void RgQuad9G9::project_to_nodes(double* ai, double* ao) const
{
    std::vector<double> b(NELN);
    for (int i = 0; i < NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j < NINT; ++j) b[i] += m_H(j, i) * ai[j];
    }

    for (int i = 0; i < NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j < NELN; ++j) ao[i] += m_Ai(i, j) * b[j];
    }
}

// Implementation for RgQuad9NI
RgQuad9NI::RgQuad9NI() : RgQuad9_(NINT, FE_QUAD9NI)
{
    double w = 1. / 9.;
    gaussPoints[0] = RgGaussPoint(-1, -1, w);
    gaussPoints[1] = RgGaussPoint( 1, -1, w);
    gaussPoints[2] = RgGaussPoint( 1,  1, w);
    gaussPoints[3] = RgGaussPoint(-1,  1, w);
    gaussPoints[4] = RgGaussPoint( 0, -1, 4 * w);
    gaussPoints[5] = RgGaussPoint( 1,  0, 4 * w);
    gaussPoints[6] = RgGaussPoint( 0,  1, 4 * w);
    gaussPoints[7] = RgGaussPoint(-1,  0, 4 * w);
    gaussPoints[8] = RgGaussPoint( 0,  0, 16 * w);
    init();
}

void RgQuad9NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
    ao[4] = ai[4];
    ao[5] = ai[5];
    ao[6] = ai[6];
    ao[7] = ai[7];
    ao[8] = ai[8];
}

