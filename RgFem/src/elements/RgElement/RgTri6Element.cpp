#include "RgTri6Element.h"
#include <cmath>

namespace RgFem {

RgTri6Element::RgTri6Element()
{
}

RgTri6Element::RgTri6Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgTri6Element::RgTri6Element(const RgTri6Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgTri6Element& RgTri6Element::operator=(const RgTri6Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgTri6Element::~RgTri6Element()
{
}

ElementType RgTri6Element::elementType() const
{
    return ElementType::FE_TRI6;
}

ElementShape RgTri6Element::elementShape() const
{
    return ElementShape::TRIANGLE;
}

ElementCategory RgTri6Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgTri6Element::getNumberOfGaussPoints() const
{
    return 6;  // 6-point Gauss quadrature for quadratic triangle
}

void RgTri6Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgTri6Element::N_corner(int cornerIdx, double r, double s) const
{
    // Quadratic shape functions for corner nodes using barycentric coordinates
    // Node 0: L1(2*L1 - 1), Node 1: L2(2*L2 - 1), Node 2: L3(2*L3 - 1)
    // Where L1 = 1-r-s, L2 = r, L3 = s
    
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (cornerIdx) {
        case 0:  // Node 0: L1(2*L1 - 1)
            return L1_val * (2.0 * L1_val - 1.0);
        case 1:  // Node 1: L2(2*L2 - 1)
            return L2_val * (2.0 * L2_val - 1.0);
        case 2:  // Node 2: L3(2*L3 - 1)
            return L3_val * (2.0 * L3_val - 1.0);
        default:
            return 0.0;
    }
}

double RgTri6Element::N_midedge(int edgeIdx, double r, double s) const
{
    // Quadratic shape functions for mid-edge nodes
    // Node 3: 4*L1*L2 (edge 0-1), Node 4: 4*L2*L3 (edge 1-2), Node 5: 4*L3*L1 (edge 2-0)
    
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (edgeIdx) {
        case 3:  // Node 3: 4*L1*L2
            return 4.0 * L1_val * L2_val;
        case 4:  // Node 4: 4*L2*L3
            return 4.0 * L2_val * L3_val;
        case 5:  // Node 5: 4*L3*L1
            return 4.0 * L3_val * L1_val;
        default:
            return 0.0;
    }
}

double RgTri6Element::dN_corner_dr(int cornerIdx, double r, double s) const
{
    // Derivatives of corner shape functions with respect to r
    // L1 = 1-r-s, so ∂L1/∂r = -1
    // L2 = r,     so ∂L2/∂r = 1
    // L3 = s,     so ∂L3/∂r = 0
    
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (cornerIdx) {
        case 0:  // ∂N0/∂r = ∂L1/∂r * (2*L1 - 1) + L1 * 2 * ∂L1/∂r = -1 * (2*L1 - 1) + L1 * 2 * (-1)
            return -(2.0 * L1_val - 1.0) - 2.0 * L1_val;
        case 1:  // ∂N1/∂r = 1 * (2*L2 - 1) + L2 * 2 * 1
            return (2.0 * L2_val - 1.0) + 2.0 * L2_val;
        case 2:  // ∂N2/∂r = 0
            return 0.0;
        default:
            return 0.0;
    }
}

double RgTri6Element::dN_corner_ds(int cornerIdx, double r, double s) const
{
    // Derivatives of corner shape functions with respect to s
    // L1 = 1-r-s, so ∂L1/∂s = -1
    // L2 = r,     so ∂L2/∂s = 0
    // L3 = s,     so ∂L3/∂s = 1
    
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (cornerIdx) {
        case 0:  // ∂N0/∂s = -1 * (2*L1 - 1) - 2*L1
            return -(2.0 * L1_val - 1.0) - 2.0 * L1_val;
        case 1:  // ∂N1/∂s = 0
            return 0.0;
        case 2:  // ∂N2/∂s = 1 * (2*L3 - 1) + L3 * 2 * 1
            return (2.0 * L3_val - 1.0) + 2.0 * L3_val;
        default:
            return 0.0;
    }
}

double RgTri6Element::dN_midedge_dr(int edgeIdx, double r, double s) const
{
    // Derivatives of mid-edge shape functions with respect to r
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (edgeIdx) {
        case 3:  // ∂(4*L1*L2)/∂r = 4 * (-1*L2 + L1*1) = 4*(L1 - L2)
            return 4.0 * (L1_val - L2_val);
        case 4:  // ∂(4*L2*L3)/∂r = 4 * (1*L3 + L2*0) = 4*L3
            return 4.0 * L3_val;
        case 5:  // ∂(4*L3*L1)/∂r = 4 * (0*L1 + L3*(-1)) = -4*L3
            return -4.0 * L3_val;
        default:
            return 0.0;
    }
}

double RgTri6Element::dN_midedge_ds(int edgeIdx, double r, double s) const
{
    // Derivatives of mid-edge shape functions with respect to s
    double L1_val = 1.0 - r - s;
    double L2_val = r;
    double L3_val = s;
    
    switch (edgeIdx) {
        case 3:  // ∂(4*L1*L2)/∂s = 4 * (-1*L2 + L1*0) = -4*L2
            return -4.0 * L2_val;
        case 4:  // ∂(4*L2*L3)/∂s = 4 * (0*L3 + L2*1) = 4*L2
            return 4.0 * L2_val;
        case 5:  // ∂(4*L3*L1)/∂s = 4 * (1*L1 + L3*(-1)) = 4*(L1 - L3)
            return 4.0 * (L1_val - L3_val);
        default:
            return 0.0;
    }
}

double RgTri6Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    if (nodeId < 3) {
        return N_corner(nodeId, r, s);
    } else {
        return N_midedge(nodeId, r, s);
    }
}

void RgTri6Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                     double& dNdr_out, double& dNds_out, double& dNdt) const
{
    dNdt = 0.0;  // No variation in z-direction for 2D element
    
    if (nodeId < 3) {
        dNdr_out = dN_corner_dr(nodeId, r, s);
        dNds_out = dN_corner_ds(nodeId, r, s);
    } else {
        dNdr_out = dN_midedge_dr(nodeId, r, s);
        dNds_out = dN_midedge_ds(nodeId, r, s);
    }
}

void RgTri6Element::evaluateCoordinates(double r, double s, double t,
                                        std::array<double, 3>& coord) const
{
    coord[0] = coord[1] = coord[2] = 0.0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        double N_val = shapeFunction(i, r, s, t);
        const auto& nodeCoord = getNodeCoordinate(i);
        coord[0] += N_val * nodeCoord[0];
        coord[1] += N_val * nodeCoord[1];
        coord[2] += N_val * nodeCoord[2];
    }
}

void RgTri6Element::evaluateJacobian(double r, double s, double t,
                                     std::array<std::array<double, 3>, 3>& J) const
{
    // Initialize Jacobian
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            J[i][j] = 0.0;
        }
    }

    // J[i][j] = ∂x_i/∂ξ_j
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr_val, dNds_val, dNdt_val;
        shapeDerivatives(node, r, s, t, dNdr_val, dNds_val, dNdt_val);
        
        const auto& coord = getNodeCoordinate(node);
        
        // Columns: ∂/∂r, ∂/∂s
        J[0][0] += dNdr_val * coord[0];
        J[1][0] += dNdr_val * coord[1];
        J[2][0] += dNdr_val * coord[2];
        
        J[0][1] += dNds_val * coord[0];
        J[1][1] += dNds_val * coord[1];
        J[2][1] += dNds_val * coord[2];
    }
}

double RgTri6Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D element: det(J) = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

double RgTri6Element::getElementArea() const
{
    // Compute area using 6-point Gauss quadrature for triangle
    // 6-point formula (Hammer, Marlowe, Stroud)
    const double gp1 = 0.091576213509770782287;
    const double gp2 = 0.816496580927725232062;
    const double gp3 = 0.108103018168070565458;
    
    // Gauss points in (r, s) coordinates
    const double r_gp[] = { gp1, gp2, gp1, gp2, gp1, gp2 };
    const double s_gp[] = { gp1, gp1, gp2, gp2, gp3, gp3 };
    const double wt = 1.0 / 6.0;
    
    double area = 0.0;
    for (int i = 0; i < 6; ++i) {
        double det = evaluateJacobianDeterminant(r_gp[i], s_gp[i], 0.0);
        area += wt * det;
    }
    
    return std::abs(area);
}

void RgTri6Element::evaluateJacobianInverse(double r, double s, double t,
                                            std::array<std::array<double, 3>, 3>& Jinv) const
{
    // For 2D element in parametric domain
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    double det = evaluateJacobianDeterminant(r, s, t);
    if (std::abs(det) < 1e-15) {
        // Singular Jacobian - set to identity
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Jinv[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return;
    }
    
    // Inverse of 2D Jacobian (embedded in 3D)
    Jinv[0][0] =  J[1][1] / det;
    Jinv[0][1] = -J[0][1] / det;
    Jinv[1][0] = -J[1][0] / det;
    Jinv[1][1] =  J[0][0] / det;
    Jinv[2][2] = 1.0;  // For out-of-plane direction
}

void RgTri6Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV (using 6-point Gauss quadrature)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 6-point Gauss quadrature
}

void RgTri6Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgTri6Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
