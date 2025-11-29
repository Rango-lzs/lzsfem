#include "RgPyramid5Element.h"
#include <cmath>

namespace RgFem {

RgPyramid5Element::RgPyramid5Element()
{
}

RgPyramid5Element::RgPyramid5Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgPyramid5Element::RgPyramid5Element(const RgPyramid5Element& other)
    : RgLinearSolid3dElement(other)
{
}

RgPyramid5Element& RgPyramid5Element::operator=(const RgPyramid5Element& other)
{
    if (this != &other) {
        RgLinearSolid3dElement::operator=(other);
    }
    return *this;
}

RgPyramid5Element::~RgPyramid5Element()
{
}

ElementType RgPyramid5Element::elementType() const
{
    return ElementType::FE_PYRAMID5;
}

ElementShape RgPyramid5Element::elementShape() const
{
    return ElementShape::PYRAMID;
}

ElementCategory RgPyramid5Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgPyramid5Element::getNumberOfGaussPoints() const
{
    return 8;  // Reduced integration for pyramid
}

void RgPyramid5Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgPyramid5Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Pyramid shape functions (isoparametric)
    // Natural coordinates: r,s in [-1,1], t in [-1,1] (t=-1 at base, t=1 at apex)
    // Base nodes 0-3 at t=-1, Apex node 4 at t=1
    
    double z = 0.5 * (1.0 - t);  // Height factor: 0 at apex, 1 at base
    
    if (nodeId < 4) {
        // Base quad nodes (bilinear in r,s)
        // Node 0: r=-1, s=-1
        // Node 1: r=1,  s=-1
        // Node 2: r=1,  s=1
        // Node 3: r=-1, s=1
        
        double N_base = 0.25 * (1.0 + r * (nodeId % 2 == 1 ? 1.0 : -1.0)) 
                              * (1.0 + s * (nodeId / 2 == 1 ? 1.0 : -1.0));
        return z * N_base;
    } else {
        // Apex node
        return 0.5 * (1.0 + t);
    }
}

void RgPyramid5Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                         double& dNdr, double& dNds, double& dNdt) const
{
    double z = 0.5 * (1.0 - t);
    
    if (nodeId < 4) {
        // Base nodes
        double r_sign = (nodeId % 2 == 1) ? 1.0 : -1.0;
        double s_sign = (nodeId / 2 == 1) ? 1.0 : -1.0;
        
        double N_base = 0.25 * (1.0 + r * r_sign) * (1.0 + s * s_sign);
        
        dNdr = 0.25 * r_sign * (1.0 + s * s_sign) * z;
        dNds = 0.25 * (1.0 + r * r_sign) * s_sign * z;
        dNdt = -0.5 * N_base;
    } else {
        // Apex node
        dNdr = 0.0;
        dNds = 0.0;
        dNdt = 0.5;
    }
}

double RgPyramid5Element::computeHeightFactor(double t) const
{
    return 0.5 * (1.0 - t);
}

void RgPyramid5Element::evaluateCoordinates(double r, double s, double t,
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

void RgPyramid5Element::evaluateJacobian(double r, double s, double t,
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
        
        J[0][2] += dNdt * coord[0];
        J[1][2] += dNdt * coord[1];
        J[2][2] += dNdt * coord[2];
    }
}

double RgPyramid5Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
         - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
         + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

void RgPyramid5Element::evaluateJacobianInverse(double r, double s, double t,
                                               std::array<std::array<double, 3>, 3>& Jinv) const
{
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
    
    Jinv[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
    Jinv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) / det;
    Jinv[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
    
    Jinv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) / det;
    Jinv[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
    Jinv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) / det;
    
    Jinv[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
    Jinv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) / det;
    Jinv[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
}

void RgPyramid5Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV using 8-point Gauss quadrature
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
}

void RgPyramid5Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgPyramid5Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
