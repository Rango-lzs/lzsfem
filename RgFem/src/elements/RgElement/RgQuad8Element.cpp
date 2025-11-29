#include "RgQuad8Element.h"
#include <cmath>

namespace RgFem {

RgQuad8Element::RgQuad8Element()
{
}

RgQuad8Element::RgQuad8Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgQuad8Element::RgQuad8Element(const RgQuad8Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgQuad8Element& RgQuad8Element::operator=(const RgQuad8Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgQuad8Element::~RgQuad8Element()
{
}

ElementType RgQuad8Element::elementType() const
{
    return ElementType::FE_QUAD8;
}

ElementShape RgQuad8Element::elementShape() const
{
    return ElementShape::QUADRILATERAL;
}

ElementCategory RgQuad8Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgQuad8Element::getNumberOfGaussPoints() const
{
    return 9;  // 3×3 Gauss quadrature for quadratic elements
}

void RgQuad8Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgQuad8Element::N_corner(int cornerIdx, double r, double s) const
{
    // Corner nodes (0,1,2,3) for serendipity quadratic element
    // Node 0: bottom-left, Node 1: bottom-right, Node 2: top-right, Node 3: top-left
    double rs = r * s;
    double r2 = r * r;
    double s2 = s * s;
    
    switch (cornerIdx) {
        case 0:  // (r=-1, s=-1): 0.25 * (1-r) * (1-s) * (r+s+1)
            return 0.25 * (1.0 - r) * (1.0 - s) * (r + s + 1.0);
        case 1:  // (r=+1, s=-1): 0.25 * (1+r) * (1-s) * (-r+s+1)
            return 0.25 * (1.0 + r) * (1.0 - s) * (-r + s + 1.0);
        case 2:  // (r=+1, s=+1): 0.25 * (1+r) * (1+s) * (r+s-1)
            return 0.25 * (1.0 + r) * (1.0 + s) * (r + s - 1.0);
        case 3:  // (r=-1, s=+1): 0.25 * (1-r) * (1+s) * (-r-s+1)
            return 0.25 * (1.0 - r) * (1.0 + s) * (-r - s + 1.0);
        default:
            return 0.0;
    }
}

double RgQuad8Element::N_midedge(int edgeIdx, double r, double s) const
{
    // Mid-edge nodes (4,5,6,7) for serendipity quadratic element
    // Node 4: bottom (s=-1), Node 5: right (r=+1), Node 6: top (s=+1), Node 7: left (r=-1)
    double r2 = r * r;
    double s2 = s * s;
    
    switch (edgeIdx) {
        case 4:  // Bottom edge (s=-1): 0.5 * (1-r²) * (1-s)
            return 0.5 * (1.0 - r2) * (1.0 - s);
        case 5:  // Right edge (r=+1): 0.5 * (1+r) * (1-s²)
            return 0.5 * (1.0 + r) * (1.0 - s2);
        case 6:  // Top edge (s=+1): 0.5 * (1-r²) * (1+s)
            return 0.5 * (1.0 - r2) * (1.0 + s);
        case 7:  // Left edge (r=-1): 0.5 * (1-r) * (1-s²)
            return 0.5 * (1.0 - r) * (1.0 - s2);
        default:
            return 0.0;
    }
}

double RgQuad8Element::dN_corner_dr(int cornerIdx, double r, double s) const
{
    // Derivatives of corner shape functions with respect to r
    double s2 = s * s;
    
    switch (cornerIdx) {
        case 0:  // ∂N0/∂r
            return 0.25 * (-(1.0 - s) * (r + s + 1.0) + (1.0 - r) * (1.0 - s));
        case 1:  // ∂N1/∂r
            return 0.25 * ((1.0 - s) * (-r + s + 1.0) + (1.0 + r) * (1.0 - s) * (-1.0));
        case 2:  // ∂N2/∂r
            return 0.25 * ((1.0 + s) * (r + s - 1.0) + (1.0 + r) * (1.0 + s));
        case 3:  // ∂N3/∂r
            return 0.25 * (-(1.0 + s) * (-r - s + 1.0) + (1.0 - r) * (1.0 + s) * (-1.0));
        default:
            return 0.0;
    }
}

double RgQuad8Element::dN_corner_ds(int cornerIdx, double r, double s) const
{
    // Derivatives of corner shape functions with respect to s
    double r2 = r * r;
    
    switch (cornerIdx) {
        case 0:  // ∂N0/∂s
            return 0.25 * (-(1.0 - r) * (r + s + 1.0) + (1.0 - r) * (1.0 - s));
        case 1:  // ∂N1/∂s
            return 0.25 * (-(1.0 + r) * (-r + s + 1.0) + (1.0 + r) * (1.0 - s));
        case 2:  // ∂N2/∂s
            return 0.25 * ((1.0 + r) * (r + s - 1.0) + (1.0 + r) * (1.0 + s));
        case 3:  // ∂N3/∂s
            return 0.25 * (-(1.0 - r) * (-r - s + 1.0) + (1.0 - r) * (1.0 + s));
        default:
            return 0.0;
    }
}

double RgQuad8Element::dN_midedge_dr(int edgeIdx, double r, double s) const
{
    // Derivatives of mid-edge shape functions with respect to r
    switch (edgeIdx) {
        case 4:  // ∂N4/∂r
            return -r * (1.0 - s);
        case 5:  // ∂N5/∂r
            return 0.5 * (1.0 - s * s);
        case 6:  // ∂N6/∂r
            return -r * (1.0 + s);
        case 7:  // ∂N7/∂r
            return -0.5 * (1.0 - s * s);
        default:
            return 0.0;
    }
}

double RgQuad8Element::dN_midedge_ds(int edgeIdx, double r, double s) const
{
    // Derivatives of mid-edge shape functions with respect to s
    switch (edgeIdx) {
        case 4:  // ∂N4/∂s
            return -0.5 * (1.0 - r * r);
        case 5:  // ∂N5/∂s
            return -s * (1.0 + r);
        case 6:  // ∂N6/∂s
            return 0.5 * (1.0 - r * r);
        case 7:  // ∂N7/∂s
            return -s * (1.0 - r);
        default:
            return 0.0;
    }
}

double RgQuad8Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    if (nodeId < 4) {
        return N_corner(nodeId, r, s);
    } else {
        return N_midedge(nodeId, r, s);
    }
}

void RgQuad8Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                      double& dNdr_out, double& dNds_out, double& dNdt) const
{
    dNdt = 0.0;  // No variation in z-direction for 2D element
    
    if (nodeId < 4) {
        dNdr_out = dN_corner_dr(nodeId, r, s);
        dNds_out = dN_corner_ds(nodeId, r, s);
    } else {
        dNdr_out = dN_midedge_dr(nodeId, r, s);
        dNds_out = dN_midedge_ds(nodeId, r, s);
    }
}

void RgQuad8Element::evaluateCoordinates(double r, double s, double t,
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

void RgQuad8Element::evaluateJacobian(double r, double s, double t,
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

double RgQuad8Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D element: det(J) = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

double RgQuad8Element::getElementArea() const
{
    // Compute area using 3×3 Gauss quadrature
    // Gauss points and weights for 3-point rule
    const double gp[] = { -0.774596669241483, 0.0, 0.774596669241483 };
    const double wt[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
    
    double area = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double det = evaluateJacobianDeterminant(gp[i], gp[j], 0.0);
            area += wt[i] * wt[j] * det;
        }
    }
    
    return std::abs(area);
}

void RgQuad8Element::evaluateJacobianInverse(double r, double s, double t,
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

void RgQuad8Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV (using 3×3 Gauss quadrature)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 3×3 Gauss points
}

void RgQuad8Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgQuad8Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
