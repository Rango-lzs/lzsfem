#include "RgTet10Element.h"
#include "elements/ElementShape/RgElementShape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <array>

namespace RgFem {

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgTet10Element::RgTet10Element()
    : RgSolid3dElement()
{
    initTraits();
}

RgTet10Element::RgTet10Element(const std::array<int, kNodeCount>& nodeIds)
    : RgSolid3dElement()
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
    initTraits();
}

RgTet10Element::RgTet10Element(const RgTet10Element& other)
    : RgSolid3dElement(other)
{
    m_gaussR = other.m_gaussR;
    m_gaussS = other.m_gaussS;
    m_gaussT = other.m_gaussT;
    m_gaussW = other.m_gaussW;
    m_jacobianInverse = other.m_jacobianInverse;
    m_jacobianDet = other.m_jacobianDet;
}

RgTet10Element::~RgTet10Element()
{
}

RgTet10Element& RgTet10Element::operator=(const RgTet10Element& other)
{
    if (this != &other) {
        RgSolid3dElement::operator=(other);
        m_gaussR = other.m_gaussR;
        m_gaussS = other.m_gaussS;
        m_gaussT = other.m_gaussT;
        m_gaussW = other.m_gaussW;
        m_jacobianInverse = other.m_jacobianInverse;
        m_jacobianDet = other.m_jacobianDet;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgTet10Element::elementType() const
{
    return ElementType::FE_TET10G4;
}

ElementShape RgTet10Element::elementShape() const
{
    return ElementShape::ET_TET10;
}

ElementCategory RgTet10Element::elementClass() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

// ============================================================================
// Element Properties
// ============================================================================

int RgTet10Element::getNumberOfGaussPoints() const
{
    return m_pTraits ? m_pTraits->m_nint : 4;
}

// ============================================================================
// Initialize Element Traits
// ============================================================================

void RgTet10Element::initTraits()
{
    RgSolidElementTraits* traits = new RgSolidElementTraits(4, kNodeCount, ET_TET10, FE_TET10G4);
    traits->init();
    
    if (m_pTraits == nullptr) {
        m_pTraits = traits;
        m_node.resize(kNodeCount);
        m_loc_node.resize(kNodeCount);
    }
    
    initializeGaussPoints();
}

void RgTet10Element::initializeGaussPoints()
{
    // 4 Gauss points for tetrahedral quadrature (full integration for quadratic elements)
    // Using standard tetrahedral quadrature points
    m_gaussR.resize(4);
    m_gaussS.resize(4);
    m_gaussT.resize(4);
    m_gaussW.resize(4);
    
    // Gauss points and weights for tetrahedral quadrature
    // Points in barycentric coordinates (L1, L2, L3, L4)
    const double a = 0.585410196627215;
    const double b = 0.138196601125011;
    
    // Gauss point 1: (a, b, b)
    m_gaussR[0] = a;
    m_gaussS[0] = b;
    m_gaussT[0] = b;
    m_gaussW[0] = 0.25;
    
    // Gauss point 2: (b, a, b)
    m_gaussR[1] = b;
    m_gaussS[1] = a;
    m_gaussT[1] = b;
    m_gaussW[1] = 0.25;
    
    // Gauss point 3: (b, b, a)
    m_gaussR[2] = b;
    m_gaussS[2] = b;
    m_gaussT[2] = a;
    m_gaussW[2] = 0.25;
    
    // Gauss point 4: (b, b, b)
    m_gaussR[3] = b;
    m_gaussS[3] = b;
    m_gaussT[3] = b;
    m_gaussW[3] = 0.25;
    
    m_jacobianInverse.resize(4);
    m_jacobianDet.resize(4);
}

// ============================================================================
// Shape Function Evaluations
// ============================================================================

void RgTet10Element::evaluateShapeFunctions(double r, double s, double t, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    
    // Quadratic tetrahedral shape functions using barycentric coordinates
    // L1 = 1 - r - s - t
    // L2 = r, L3 = s, L4 = t
    
    double L1 = 1.0 - r - s - t;
    double L2 = r;
    double L3 = s;
    double L4 = t;
    
    // Corner nodes (quadratic)
    N[0] = L1 * (2.0 * L1 - 1.0);  // Node 0: (0,0,0)
    N[1] = L2 * (2.0 * L2 - 1.0);  // Node 1: (1,0,0)
    N[2] = L3 * (2.0 * L3 - 1.0);  // Node 2: (0,1,0)
    N[3] = L4 * (2.0 * L4 - 1.0);  // Node 3: (0,0,1)
    
    // Mid-edge nodes (linear products)
    N[4] = 4.0 * L1 * L2;  // Edge 0-1
    N[5] = 4.0 * L2 * L3;  // Edge 1-2
    N[6] = 4.0 * L3 * L1;  // Edge 2-0
    N[7] = 4.0 * L1 * L4;  // Edge 0-3
    N[8] = 4.0 * L2 * L4;  // Edge 1-3
    N[9] = 4.0 * L3 * L4;  // Edge 2-3
}

void RgTet10Element::evaluateShapeDerivatives(double r, double s, double t,
                                              std::vector<double>& dN_dr,
                                              std::vector<double>& dN_ds,
                                              std::vector<double>& dN_dt) const
{
    dN_dr.assign(kNodeCount, 0.0);
    dN_ds.assign(kNodeCount, 0.0);
    dN_dt.assign(kNodeCount, 0.0);
    
    double L1 = 1.0 - r - s - t;
    double L2 = r;
    double L3 = s;
    double L4 = t;
    
    // Derivatives of L_i
    // dL1/dr = -1, dL1/ds = -1, dL1/dt = -1
    // dL2/dr = 1,  dL2/ds = 0,  dL2/dt = 0
    // dL3/dr = 0,  dL3/ds = 1,  dL3/dt = 0
    // dL4/dr = 0,  dL4/ds = 0,  dL4/dt = 1
    
    // Corner nodes
    // dN0 = d(L1(2L1-1)) = 2L1(dL1) + (2L1-1)(dL1) = (4L1-1)(dL1)
    dN_dr[0] = (4.0 * L1 - 1.0) * (-1.0);
    dN_ds[0] = (4.0 * L1 - 1.0) * (-1.0);
    dN_dt[0] = (4.0 * L1 - 1.0) * (-1.0);
    
    // dN1 = d(L2(2L2-1)) = (4L2-1)(dL2)
    dN_dr[1] = (4.0 * L2 - 1.0) * 1.0;
    dN_ds[1] = 0.0;
    dN_dt[1] = 0.0;
    
    // dN2 = d(L3(2L3-1)) = (4L3-1)(dL3)
    dN_dr[2] = 0.0;
    dN_ds[2] = (4.0 * L3 - 1.0) * 1.0;
    dN_dt[2] = 0.0;
    
    // dN3 = d(L4(2L4-1)) = (4L4-1)(dL4)
    dN_dr[3] = 0.0;
    dN_ds[3] = 0.0;
    dN_dt[3] = (4.0 * L4 - 1.0) * 1.0;
    
    // Mid-edge nodes: d(4*L_i*L_j) = 4(L_i*dL_j + L_j*dL_i)
    
    // N4 = 4*L1*L2
    dN_dr[4] = 4.0 * (L1 * 1.0 + L2 * (-1.0));  // 4(L1 - L2)
    dN_ds[4] = 4.0 * (L1 * 0.0 + L2 * (-1.0));  // -4*L2
    dN_dt[4] = 4.0 * (L1 * 0.0 + L2 * (-1.0));  // -4*L2
    
    // N5 = 4*L2*L3
    dN_dr[5] = 4.0 * (L2 * 0.0 + L3 * 1.0);    // 4*L3
    dN_ds[5] = 4.0 * (L2 * 1.0 + L3 * 0.0);    // 4*L2
    dN_dt[5] = 0.0;
    
    // N6 = 4*L3*L1
    dN_dr[6] = 4.0 * (L3 * (-1.0) + L1 * 0.0); // -4*L3
    dN_ds[6] = 4.0 * (L3 * 0.0 + L1 * 1.0);    // 4*L1
    dN_dt[6] = 4.0 * (L3 * (-1.0) + L1 * 0.0); // -4*L3
    
    // N7 = 4*L1*L4
    dN_dr[7] = 4.0 * (L1 * 0.0 + L4 * (-1.0)); // -4*L4
    dN_ds[7] = 4.0 * (L1 * 0.0 + L4 * (-1.0)); // -4*L4
    dN_dt[7] = 4.0 * (L1 * 1.0 + L4 * (-1.0)); // 4(L1 - L4)
    
    // N8 = 4*L2*L4
    dN_dr[8] = 4.0 * (L2 * 0.0 + L4 * 1.0);    // 4*L4
    dN_ds[8] = 0.0;
    dN_dt[8] = 4.0 * (L2 * 1.0 + L4 * 0.0);    // 4*L2
    
    // N9 = 4*L3*L4
    dN_dr[9] = 0.0;
    dN_ds[9] = 4.0 * (L3 * 0.0 + L4 * 1.0);    // 4*L4
    dN_dt[9] = 4.0 * (L3 * 1.0 + L4 * 0.0);    // 4*L3
}

void RgTet10Element::evaluateShapeDerivatives2(double r, double s, double t,
                                               std::vector<double>& d2N_drr,
                                               std::vector<double>& d2N_dss,
                                               std::vector<double>& d2N_dtt,
                                               std::vector<double>& d2N_drs,
                                               std::vector<double>& d2N_dst,
                                               std::vector<double>& d2N_drt) const
{
    d2N_drr.assign(kNodeCount, 0.0);
    d2N_dss.assign(kNodeCount, 0.0);
    d2N_dtt.assign(kNodeCount, 0.0);
    d2N_drs.assign(kNodeCount, 0.0);
    d2N_dst.assign(kNodeCount, 0.0);
    d2N_drt.assign(kNodeCount, 0.0);
    
    // For quadratic elements, second derivatives are constant
    // d2(N)/dr2 for corner nodes: d(4Li-1)(dLi/dr) = 4(dLi/dr)(dLi/dr)
    
    double L1 = 1.0 - r - s - t;
    
    // Corner nodes (4Li - 1)
    // d2N0/dr2 = 4*(-1)*(-1) = 4
    d2N_drr[0] = 4.0;
    d2N_dss[0] = 4.0;
    d2N_dtt[0] = 4.0;
    d2N_drs[0] = 4.0;
    d2N_dst[0] = 4.0;
    d2N_drt[0] = 4.0;
    
    // d2N1/dr2 = 4*1*1 = 4
    d2N_drr[1] = 4.0;
    d2N_dss[1] = 0.0;
    d2N_dtt[1] = 0.0;
    d2N_drs[1] = 0.0;
    d2N_dst[1] = 0.0;
    d2N_drt[1] = 0.0;
    
    // d2N2/ds2 = 4*1*1 = 4
    d2N_drr[2] = 0.0;
    d2N_dss[2] = 4.0;
    d2N_dtt[2] = 0.0;
    d2N_drs[2] = 0.0;
    d2N_dst[2] = 0.0;
    d2N_drt[2] = 0.0;
    
    // d2N3/dt2 = 4*1*1 = 4
    d2N_drr[3] = 0.0;
    d2N_dss[3] = 0.0;
    d2N_dtt[3] = 4.0;
    d2N_drs[3] = 0.0;
    d2N_dst[3] = 0.0;
    d2N_drt[3] = 0.0;
    
    // Mid-edge nodes: 4(Li*dLj + Lj*dLi)
    // Second derivatives: 4(dLi/dx*dLj/dx + dLi/dx*dLj/dx) = 4*2(dLi/dx)(dLj/dx)
    
    // N4 = 4(L1*L2): d2/dr2 = 4*2*(-1)*1 = -8
    d2N_drr[4] = -8.0;
    d2N_dss[4] = 0.0;
    d2N_dtt[4] = 0.0;
    d2N_drs[4] = 0.0;
    d2N_dst[4] = 0.0;
    d2N_drt[4] = 0.0;
    
    // N5 = 4(L2*L3): d2/drs = 4*2*1*1 = 8
    d2N_drr[5] = 0.0;
    d2N_dss[5] = 0.0;
    d2N_dtt[5] = 0.0;
    d2N_drs[5] = 8.0;
    d2N_dst[5] = 0.0;
    d2N_drt[5] = 0.0;
    
    // N6 = 4(L3*L1): d2/ds2 = 4*2*1*(-1) = -8, d2/drs = 4*2*1*(-1) = -8
    d2N_drr[6] = 0.0;
    d2N_dss[6] = -8.0;
    d2N_dtt[6] = 0.0;
    d2N_drs[6] = -8.0;
    d2N_dst[6] = 0.0;
    d2N_drt[6] = 0.0;
    
    // N7 = 4(L1*L4): d2/dt2 = 4*2*(-1)*1 = -8, d2/drt = 4*2*(-1)*1 = -8
    d2N_drr[7] = 0.0;
    d2N_dss[7] = 0.0;
    d2N_dtt[7] = -8.0;
    d2N_drs[7] = 0.0;
    d2N_dst[7] = -8.0;
    d2N_drt[7] = -8.0;
    
    // N8 = 4(L2*L4): d2/drt = 4*2*1*1 = 8
    d2N_drr[8] = 0.0;
    d2N_dss[8] = 0.0;
    d2N_dtt[8] = 0.0;
    d2N_drs[8] = 0.0;
    d2N_dst[8] = 0.0;
    d2N_drt[8] = 8.0;
    
    // N9 = 4(L3*L4): d2/dst = 4*2*1*1 = 8
    d2N_drr[9] = 0.0;
    d2N_dss[9] = 0.0;
    d2N_dtt[9] = 0.0;
    d2N_drs[9] = 0.0;
    d2N_dst[9] = 8.0;
    d2N_drt[9] = 0.0;
}

// ============================================================================
// Geometric Operations
// ============================================================================

Vector3d RgTet10Element::evaluateCoordinates(const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d coords(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            coords += N[i] * node->getPosition();
        }
    }
    return coords;
}

Matrix3d RgTet10Element::evaluateJacobian(const Vector3d& naturalCoord) const
{
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d J(0, 0, 0, 0, 0, 0, 0, 0, 0);
    
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            Vector3d pos = node->getPosition();
            J.m[0][0] += dN_dr[i] * pos.x;
            J.m[0][1] += dN_dr[i] * pos.y;
            J.m[0][2] += dN_dr[i] * pos.z;
            J.m[1][0] += dN_ds[i] * pos.x;
            J.m[1][1] += dN_ds[i] * pos.y;
            J.m[1][2] += dN_ds[i] * pos.z;
            J.m[2][0] += dN_dt[i] * pos.x;
            J.m[2][1] += dN_dt[i] * pos.y;
            J.m[2][2] += dN_dt[i] * pos.z;
        }
    }
    return J;
}

double RgTet10Element::evaluateJacobianDeterminant(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.det();
}

Matrix3d RgTet10Element::evaluateJacobianInverse(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.inverse();
}

// ============================================================================
// Physical Field Evaluations
// ============================================================================

Vector3d RgTet10Element::evaluateField(const Vector3d* nodeValues, const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d field(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        field += N[i] * nodeValues[i];
    }
    return field;
}

double RgTet10Element::evaluateScalarField(const double* nodeValues, const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    double field = 0.0;
    for (int i = 0; i < kNodeCount; ++i) {
        field += N[i] * nodeValues[i];
    }
    return field;
}

// ============================================================================
// B-Matrix Computation
// ============================================================================

void RgTet10Element::computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const
{
    int ndofs = kNodeCount * 3;
    B.resize(6, ndofs);
    B.zero();
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    JinvT = JinvT.transpose();
    
    // Transform shape derivatives from natural to physical coordinates
    for (int i = 0; i < kNodeCount; ++i) {
        double dN_dx = JinvT.m[0][0] * dN_dr[i] + JinvT.m[0][1] * dN_ds[i] + JinvT.m[0][2] * dN_dt[i];
        double dN_dy = JinvT.m[1][0] * dN_dr[i] + JinvT.m[1][1] * dN_ds[i] + JinvT.m[1][2] * dN_dt[i];
        double dN_dz = JinvT.m[2][0] * dN_dr[i] + JinvT.m[2][1] * dN_ds[i] + JinvT.m[2][2] * dN_dt[i];
        
        // Strain-displacement matrix for node i (3 DOFs: u, v, w)
        // Voigt notation: [εxx, εyy, εzz, γxy, γyz, γzx]
        B[0][3*i]     = dN_dx;   // εxx = ∂u/∂x
        B[1][3*i+1]   = dN_dy;   // εyy = ∂v/∂y
        B[2][3*i+2]   = dN_dz;   // εzz = ∂w/∂z
        B[3][3*i]     = dN_dy;   // γxy = ∂u/∂y
        B[3][3*i+1]   = dN_dx;   // γxy = ∂v/∂x
        B[4][3*i+1]   = dN_dz;   // γyz = ∂v/∂z
        B[4][3*i+2]   = dN_dy;   // γyz = ∂w/∂y
        B[5][3*i]     = dN_dz;   // γzx = ∂u/∂z
        B[5][3*i+2]   = dN_dx;   // γzx = ∂w/∂x
    }
}

// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgTet10Element::calculateStiffnessMatrix(Matrix& K) const
{
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Get material point properties
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    // Get material matrix (D matrix for isotropic linear elasticity)
    Matrix D;
    matPt->m_pMat->getConstitutiveMatrix(D);
    
    // Integration over 4 Gauss points
    for (int igp = 0; igp < 4; ++igp) {
        Vector3d natCoord(m_gaussR[igp], m_gaussS[igp], m_gaussT[igp]);
        
        Matrix B;
        computeBMatrix(natCoord, B);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[igp] * jdet;
        
        // K += B^T * D * B * weight
        Matrix BtD(B.rows(), D.cols());
        BtD.zero();
        
        for (int i = 0; i < B.rows(); ++i) {
            for (int j = 0; j < D.cols(); ++j) {
                for (int k = 0; k < B.cols(); ++k) {
                    BtD[i][j] += B[i][k] * D[k][j];
                }
            }
        }
        
        for (int i = 0; i < ndofs; ++i) {
            for (int j = 0; j < ndofs; ++j) {
                for (int k = 0; k < B.rows(); ++k) {
                    K[i][j] += B[k][i] * BtD[k][j] * weight;
                }
            }
        }
    }
}

void RgTet10Element::calculateMassMatrix(Matrix& M) const
{
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Get material point properties
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    double density = matPt->m_pMat->getDensity();
    
    // Consistent mass matrix: M = ρ * ∫ N^T * N dV
    for (int igp = 0; igp < 4; ++igp) {
        Vector3d natCoord(m_gaussR[igp], m_gaussS[igp], m_gaussT[igp]);
        
        std::vector<double> N;
        evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = density * m_gaussW[igp] * jdet;
        
        for (int i = 0; i < kNodeCount; ++i) {
            for (int j = 0; j < kNodeCount; ++j) {
                double Nij = N[i] * N[j] * weight;
                // Diagonal blocks for each DOF
                M[3*i][3*j]       += Nij;
                M[3*i+1][3*j+1]   += Nij;
                M[3*i+2][3*j+2]   += Nij;
            }
        }
    }
}

void RgTet10Element::calculateDampingMatrix(Matrix& C) const
{
    int ndofs = kNodeCount * 3;
    C.resize(ndofs, ndofs);
    C.zero();
    
    // Rayleigh damping: C = a0*M + a1*K
    Matrix M, K;
    calculateMassMatrix(M);
    calculateStiffnessMatrix(K);
    
    // Default Rayleigh coefficients
    double a0 = 0.0;
    double a1 = 0.0;
    
    for (int i = 0; i < ndofs; ++i) {
        for (int j = 0; j < ndofs; ++j) {
            C[i][j] = a0 * M[i][j] + a1 * K[i][j];
        }
    }
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgTet10Element::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    if (!matPt.m_pMat) return;
    
    // Compute strain from displacement
    Matrix3ds strain;
    calculateStrain(matPt, strain);
    
    // Get material constitutive relation
    matPt.m_pMat->getStress(matPt, strain, stress);
}

void RgTet10Element::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    strain.zero();
    
    // Small strain: ε = 0.5(∇u + ∇u^T)
    // Using B-matrix: strain = B * u
}

// ============================================================================
// Face and Edge Operations
// ============================================================================

void RgTet10Element::getFaceNodeIds(int faceId, std::array<int, 6>& faceNodes) const
{
    // Triangular faces with 6 nodes each (3 corners + 3 mid-edges)
    // Face numbering: 0:(0,1,2), 1:(0,1,3), 2:(0,2,3), 3:(1,2,3)
    
    switch (faceId) {
        case 0:
            faceNodes = {0, 1, 2, 4, 5, 6};  // Corner nodes 0,1,2; Mid-edges 0-1,1-2,2-0
            break;
        case 1:
            faceNodes = {0, 1, 3, 4, 8, 7};  // Corner nodes 0,1,3; Mid-edges 0-1,1-3,3-0
            break;
        case 2:
            faceNodes = {0, 2, 3, 6, 9, 7};  // Corner nodes 0,2,3; Mid-edges 0-2,2-3,3-0
            break;
        case 3:
            faceNodes = {1, 2, 3, 5, 9, 8};  // Corner nodes 1,2,3; Mid-edges 1-2,2-3,3-1
            break;
        default:
            faceNodes.fill(-1);
    }
}

void RgTet10Element::getEdgeNodeIds(int edgeId, std::array<int, 3>& edgeNodes) const
{
    // Edges with 3 nodes each (2 corners + 1 mid-point)
    // Edge numbering: 0:(0-1), 1:(1-2), 2:(2-0), 3:(0-3), 4:(1-3), 5:(2-3)
    
    switch (edgeId) {
        case 0:
            edgeNodes = {0, 1, 4};  // Edge 0-1 with mid-node 4
            break;
        case 1:
            edgeNodes = {1, 2, 5};  // Edge 1-2 with mid-node 5
            break;
        case 2:
            edgeNodes = {2, 0, 6};  // Edge 2-0 with mid-node 6
            break;
        case 3:
            edgeNodes = {0, 3, 7};  // Edge 0-3 with mid-node 7
            break;
        case 4:
            edgeNodes = {1, 3, 8};  // Edge 1-3 with mid-node 8
            break;
        case 5:
            edgeNodes = {2, 3, 9};  // Edge 2-3 with mid-node 9
            break;
        default:
            edgeNodes.fill(-1);
    }
}

// ============================================================================
// Loading Operations
// ============================================================================

void RgTet10Element::applyBodyForce(const Vector3d& force, Vector& F) const
{
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    // F = ∫ N^T * f_body dV
    for (int igp = 0; igp < 4; ++igp) {
        Vector3d natCoord(m_gaussR[igp], m_gaussS[igp], m_gaussT[igp]);
        
        std::vector<double> N;
        evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[igp] * jdet;
        
        for (int i = 0; i < kNodeCount; ++i) {
            F[3*i]   += N[i] * force.x * weight;
            F[3*i+1] += N[i] * force.y * weight;
            F[3*i+2] += N[i] * force.z * weight;
        }
    }
}

void RgTet10Element::applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const
{
    // Apply traction on triangular face
    std::array<int, 6> faceNodeIds;
    getFaceNodeIds(faceId, faceNodeIds);
    
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    // Surface traction: F = ∫_S N^T * t dS
    // For triangular elements on tet faces, use 3-point Gauss quadrature
    
    // This requires transformation to face local coordinates
    // Simplified implementation: distribute to face nodes
    for (int i = 0; i < 6; ++i) {
        int nodeId = faceNodeIds[i];
        F[3*nodeId]   += traction.x / 6.0;
        F[3*nodeId+1] += traction.y / 6.0;
        F[3*nodeId+2] += traction.z / 6.0;
    }
}

void RgTet10Element::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    if (nodeId >= 0 && nodeId < kNodeCount) {
        F[3*nodeId]   = force.x;
        F[3*nodeId+1] = force.y;
        F[3*nodeId+2] = force.z;
    }
}

// ============================================================================
// Material Point Access
// ============================================================================

FEMaterialPoint* RgTet10Element::getMaterialPoint(int gaussPtId)
{
    if (gaussPtId >= 0 && gaussPtId < m_pTraits->m_nint && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

const FEMaterialPoint* RgTet10Element::getMaterialPoint(int gaussPtId) const
{
    if (gaussPtId >= 0 && gaussPtId < m_pTraits->m_nint && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

// ============================================================================
// Serialization
// ============================================================================

void RgTet10Element::Serialize(DumpStream& ar)
{
    RgSolid3dElement::Serialize(ar);
    
    if (ar.IsSaving()) {
        ar << (int)m_gaussR.size();
        for (double v : m_gaussR) ar << v;
        for (double v : m_gaussS) ar << v;
        for (double v : m_gaussT) ar << v;
        for (double v : m_gaussW) ar << v;
    } else {
        int size = 0;
        ar >> size;
        m_gaussR.resize(size);
        m_gaussS.resize(size);
        m_gaussT.resize(size);
        m_gaussW.resize(size);
        for (int i = 0; i < size; ++i) {
            ar >> m_gaussR[i] >> m_gaussS[i] >> m_gaussT[i] >> m_gaussW[i];
        }
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

bool RgTet10Element::isValidNaturalCoordinate(const Vector3d& naturalCoord) const
{
    // Valid region: r, s, t >= 0 and r + s + t <= 1
    return (naturalCoord.x >= -1e-10 && naturalCoord.y >= -1e-10 && naturalCoord.z >= -1e-10 &&
            (naturalCoord.x + naturalCoord.y + naturalCoord.z) <= 1.0 + 1e-10);
}

double RgTet10Element::getCharacteristicLength() const
{
    // Average edge length
    double totalLength = 0.0;
    int edgeCount = 0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        for (int j = i + 1; j < kNodeCount; ++j) {
            FENode* ni = getNode(i);
            FENode* nj = getNode(j);
            if (ni && nj) {
                Vector3d diff = ni->getPosition() - nj->getPosition();
                totalLength += diff.Length();
                ++edgeCount;
            }
        }
    }
    
    return edgeCount > 0 ? totalLength / edgeCount : 1.0;
}

double RgTet10Element::getVolume() const
{
    // Tetrahedral volume from 3 edges from node 0
    FENode* n0 = getNode(0);
    FENode* n1 = getNode(1);
    FENode* n2 = getNode(2);
    FENode* n3 = getNode(3);
    
    if (!n0 || !n1 || !n2 || !n3) return 0.0;
    
    Vector3d v1 = n1->getPosition() - n0->getPosition();
    Vector3d v2 = n2->getPosition() - n0->getPosition();
    Vector3d v3 = n3->getPosition() - n0->getPosition();
    
    // Volume = |det(v1, v2, v3)| / 6
    double det = v1.x * (v2.y * v3.z - v2.z * v3.y)
               - v1.y * (v2.x * v3.z - v2.z * v3.x)
               + v1.z * (v2.x * v3.y - v2.y * v3.x);
    
    return fabs(det) / 6.0;
}

} // namespace RgFem
