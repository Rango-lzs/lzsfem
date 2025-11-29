#include "RgTri3Element.h"
#include <cmath>

namespace RgFem {

RgTri3Element::RgTri3Element()
{
}

RgTri3Element::RgTri3Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgTri3Element::RgTri3Element(const RgTri3Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgTri3Element& RgTri3Element::operator=(const RgTri3Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgTri3Element::~RgTri3Element()
{
}

ElementType RgTri3Element::elementType() const
{
    return ElementType::FE_TRI3;
}

ElementShape RgTri3Element::elementShape() const
{
    return ElementShape::TRIANGLE;
}

ElementCategory RgTri3Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgTri3Element::getNumberOfGaussPoints() const
{
    return 1;  // 1-point Gauss for constant strain
}

void RgTri3Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgTri3Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Linear triangular shape functions (barycentric coordinates)
    // Node 0: (1-r-s), Node 1: r, Node 2: s
    switch (nodeId) {
        case 0: return 1.0 - r - s;
        case 1: return r;
        case 2: return s;
        default: return 0.0;
    }
}

void RgTri3Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                     double& dNdr, double& dNds, double& dNdt) const
{
    // For 2D element, dNdt = 0 (no variation in thickness direction if treated as 3D)
    switch (nodeId) {
        case 0:
            dNdr = -1.0; dNds = -1.0; dNdt = 0.0;
            break;
        case 1:
            dNdr = 1.0; dNds = 0.0; dNdt = 0.0;
            break;
        case 2:
            dNdr = 0.0; dNds = 1.0; dNdt = 0.0;
            break;
        default:
            dNdr = 0.0; dNds = 0.0; dNdt = 0.0;
    }
}

void RgTri3Element::evaluateCoordinates(double r, double s, double t,
                                        std::array<double, 3>& coord) const
{
    coord[0] = coord[1] = coord[2] = 0.0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        double N = shapeFunction(i, r, s, t);
        const auto& nodeCoord = getNodeCoordinate(i);
        coord[0] += N * nodeCoord[0];
        coord[1] += N * nodeCoord[1];
        coord[2] += N * nodeCoord[2];
    }
}

void RgTri3Element::evaluateJacobian(double r, double s, double t,
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
        double dNdr, dNds, dNdt;
        shapeDerivatives(node, r, s, t, dNdr, dNds, dNdt);
        
        const auto& coord = getNodeCoordinate(node);
        
        J[0][0] += dNdr * coord[0];
        J[1][0] += dNdr * coord[1];
        J[2][0] += dNdr * coord[2];
        
        J[0][1] += dNds * coord[0];
        J[1][1] += dNds * coord[1];
        J[2][1] += dNds * coord[2];
    }
    
    // For 2D element in 3D space, set up in-plane Jacobian
    // J[2][2] not directly used in 2D context
}

double RgTri3Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D triangular element, det(J) is the 2D determinant
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    const auto& x2 = getNodeCoordinate(2);
    
    // Vectors from node 0 to nodes 1 and 2
    double v1_x = x1[0] - x0[0];
    double v1_y = x1[1] - x0[1];
    double v2_x = x2[0] - x0[0];
    double v2_y = x2[1] - x0[1];
    
    // Jacobian determinant = 2 * area of triangle
    return v1_x * v2_y - v1_y * v2_x;
}

double RgTri3Element::getElementArea() const
{
    return 0.5 * std::abs(evaluateJacobianDeterminant(0.0, 0.0, 0.0));
}

void RgTri3Element::evaluateJacobianInverse(double r, double s, double t,
                                           std::array<std::array<double, 3>, 3>& Jinv) const
{
    // For 2D element
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    double det = evaluateJacobianDeterminant(r, s, t);
    if (std::abs(det) < 1e-15) {
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

void RgTri3Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV (using single Gauss point for linear element)
    int ndofs = kNodeCount * 2;  // 2D element
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
}

void RgTri3Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgTri3Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
