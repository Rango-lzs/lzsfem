#include "RgQuad4Element.h"
#include <cmath>

namespace RgFem {

RgQuad4Element::RgQuad4Element()
{
}

RgQuad4Element::RgQuad4Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgQuad4Element::RgQuad4Element(const RgQuad4Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgQuad4Element& RgQuad4Element::operator=(const RgQuad4Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgQuad4Element::~RgQuad4Element()
{
}

ElementType RgQuad4Element::elementType() const
{
    return ElementType::FE_QUAD4;
}

ElementShape RgQuad4Element::elementShape() const
{
    return ElementShape::QUADRILATERAL;
}

ElementCategory RgQuad4Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgQuad4Element::getNumberOfGaussPoints() const
{
    return 4;  // 2×2 Gauss quadrature
}

void RgQuad4Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgQuad4Element::N(int nodeId, double r, double s) const
{
    // Bilinear shape functions
    switch (nodeId) {
        case 0: return 0.25 * (1.0 - r) * (1.0 - s);
        case 1: return 0.25 * (1.0 + r) * (1.0 - s);
        case 2: return 0.25 * (1.0 + r) * (1.0 + s);
        case 3: return 0.25 * (1.0 - r) * (1.0 + s);
        default: return 0.0;
    }
}

double RgQuad4Element::dNdr(int nodeId, double r, double s) const
{
    // Derivatives with respect to r
    switch (nodeId) {
        case 0: return -0.25 * (1.0 - s);
        case 1: return  0.25 * (1.0 - s);
        case 2: return  0.25 * (1.0 + s);
        case 3: return -0.25 * (1.0 + s);
        default: return 0.0;
    }
}

double RgQuad4Element::dNds(int nodeId, double r, double s) const
{
    // Derivatives with respect to s
    switch (nodeId) {
        case 0: return -0.25 * (1.0 - r);
        case 1: return -0.25 * (1.0 + r);
        case 2: return  0.25 * (1.0 + r);
        case 3: return  0.25 * (1.0 - r);
        default: return 0.0;
    }
}

double RgQuad4Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    return N(nodeId, r, s);
}

void RgQuad4Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                      double& dNdr_out, double& dNds_out, double& dNdt) const
{
    dNdr_out = dNdr(nodeId, r, s);
    dNds_out = dNds(nodeId, r, s);
    dNdt = 0.0;  // No variation in z-direction for 2D element
}

void RgQuad4Element::evaluateCoordinates(double r, double s, double t,
                                         std::array<double, 3>& coord) const
{
    coord[0] = coord[1] = coord[2] = 0.0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        double N_val = N(i, r, s);
        const auto& nodeCoord = getNodeCoordinate(i);
        coord[0] += N_val * nodeCoord[0];
        coord[1] += N_val * nodeCoord[1];
        coord[2] += N_val * nodeCoord[2];
    }
}

void RgQuad4Element::evaluateJacobian(double r, double s, double t,
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
        double dNdr_val = dNdr(node, r, s);
        double dNds_val = dNds(node, r, s);
        
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

double RgQuad4Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D element: det(J) = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

double RgQuad4Element::getElementArea() const
{
    // Compute area using 2×2 Gauss quadrature
    // Gauss points and weights
    const double gp[] = { -0.577350269189626, 0.577350269189626 };  // ±1/√3
    const double wt[] = { 1.0, 1.0 };
    
    double area = 0.0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            double det = evaluateJacobianDeterminant(gp[i], gp[j], 0.0);
            area += wt[i] * wt[j] * det;
        }
    }
    
    return std::abs(area);
}

void RgQuad4Element::evaluateJacobianInverse(double r, double s, double t,
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

void RgQuad4Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV (using 2×2 Gauss quadrature)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 2×2 Gauss points
}

void RgQuad4Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgQuad4Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
