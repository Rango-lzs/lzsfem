#include "RgSolid3dElement.h"

// return spatial dimension (3 for 3D solids)
int RgSolid3dElement::dim()
{
    return 3;
}

// 高斯点的形函数相关计算
// n : the n-th gauss point
std::vector<double> RgSolid3dElement::gaussPoint(int n) 
{
    return ((RgSolidElementTraits*)(m_pTraits))->gaussPoint(n);
}

// returns weights of integration points as a vector
std::vector<double> RgSolid3dElement::GaussWeights() const
{
    return ((RgSolidElementTraits*)(m_pTraits))->gw;
}

// dH/dr[n] -- return shape function derivative arrays as vectors
// n : the n-th gauss point
// return : N shape function derivative
std::vector<double> RgSolid3dElement::Gr(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->m_Gr[n];
}

// shape function derivative to r
std::vector<double> RgSolid3dElement::Gs(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->m_Gs[n];
}

// shape function derivative to s
std::vector<double> RgSolid3dElement::Gt(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->m_Gt[n];
}

// dH2/dr2[n] -- second derivatives as vectors
std::vector<double> RgSolid3dElement::Grr(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Grr[n];
}

std::vector<double> RgSolid3dElement::Gsr(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gsr[n];
}

std::vector<double> RgSolid3dElement::Gtr(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gtr[n];
}

std::vector<double> RgSolid3dElement::Grs(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Grs[n];
}

std::vector<double> RgSolid3dElement::Gss(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gss[n];
}

std::vector<double> RgSolid3dElement::Gts(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gts[n];
}

std::vector<double> RgSolid3dElement::Grt(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Grt[n];
}

std::vector<double> RgSolid3dElement::Gst(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gst[n];
}

std::vector<double> RgSolid3dElement::Gtt(int n) const
{
    return ((RgSolidElementTraits*)(m_pTraits))->Gtt[n];
}

//! values of shape functions (unchanged API)  这些接口计算任意点的形函数
void RgSolid3dElement::shape_fnc(double* H, double r, double s, double t) const
{
    ((RgSolidElementTraits*)(m_pTraits))->shape_fnc(H, r, s, t);
}

//! values of shape function derivatives (unchanged API)
void RgSolid3dElement::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const
{
    ((RgSolidElementTraits*)(m_pTraits))->shape_deriv(Hr, Hs, Ht, r, s, t);
}

//! values of shape function second derivatives (unchanged API)
void RgSolid3dElement::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                                     double t) const
{
    ((RgSolidElementTraits*)(m_pTraits))->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t);
}
                                    
// Additional virtual methods for 3D solid elements
int RgSolid3dElement::getNumberOfGaussPoints() const
{
    return m_pTraits ? m_pTraits->m_nint : 0;
}

void RgSolid3dElement::initTraits() 
{
    RgSolidElement::initTraits();
}

// Tangent stiffness matrix calculation
void RgSolid3dElement::calculateTangentStiffnessMatrix(Matrix& Kt) const
{
    // For linear elements, tangent stiffness equals regular stiffness
    calculateStiffnessMatrix(Kt);
}

// Geometric stiffness matrix calculation
void RgSolid3dElement::computeGeometricStiffness(Matrix& Kg) const
{
    // Default implementation: zero geometric stiffness matrix
    Kg.setZero();
}