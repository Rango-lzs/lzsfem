#include "RgWedge6Element.h"
#include <cmath>

namespace RgFem {

RgWedge6Element::RgWedge6Element()
{
}

RgWedge6Element::RgWedge6Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgWedge6Element::RgWedge6Element(const RgWedge6Element& other)
    : RgLinearSolid3dElement(other)
{
}

RgWedge6Element& RgWedge6Element::operator=(const RgWedge6Element& other)
{
    if (this != &other) {
        RgLinearSolid3dElement::operator=(other);
    }
    return *this;
}

RgWedge6Element::~RgWedge6Element()
{
}

ElementType RgWedge6Element::elementType() const
{
    return ElementType::FE_WEDGE6;
}

ElementShape RgWedge6Element::elementShape() const
{
    return ElementShape::WEDGE;
}

ElementCategory RgWedge6Element::elementClass() const
{
    return ElementCategory::SOLID;
}

int RgWedge6Element::getNumberOfGaussPoints() const
{
    return 6;  // 3-point Gauss for triangle x 2-point for z
}

void RgWedge6Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgWedge6Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Wedge shape function = Triangular(r,s) * Linear(t)
    // Node numbering (0-5):
    // Bottom (t=-1): 0(0,0,-1), 1(1,0,-1), 2(0,1,-1)
    // Top (t=1):     3(0,0,1),  4(1,0,1),  5(0,1,1)
    
    if (nodeId < 3) {
        // Bottom nodes
        double N_tri = N_tri(nodeId, r, s);
        double N_z = 0.5 * (1.0 - t);
        return N_tri * N_z;
    } else {
        // Top nodes
        double N_tri = N_tri(nodeId - 3, r, s);
        double N_z = 0.5 * (1.0 + t);
        return N_tri * N_z;
    }
}

double RgWedge6Element::N_tri(int triNodeId, double r, double s) const
{
    // Linear triangular shape functions
    // Node 0: (1-r-s), Node 1: r, Node 2: s
    switch (triNodeId) {
        case 0: return 1.0 - r - s;
        case 1: return r;
        case 2: return s;
        default: return 0.0;
    }
}

void RgWedge6Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                       double& dNdr, double& dNds, double& dNdt) const
{
    double dN_tri_r, dN_tri_s;
    int tri_id = nodeId % 3;
    dN_tri_dr(tri_id, r, s, dN_tri_r, dN_tri_s);

    if (nodeId < 3) {
        // Bottom nodes: N = N_tri * 0.5*(1-t)
        dNdr = dN_tri_r * 0.5 * (1.0 - t);
        dNds = dN_tri_s * 0.5 * (1.0 - t);
        dNdt = -0.5 * N_tri(tri_id, r, s);
    } else {
        // Top nodes: N = N_tri * 0.5*(1+t)
        dNdr = dN_tri_r * 0.5 * (1.0 + t);
        dNds = dN_tri_s * 0.5 * (1.0 + t);
        dNdt = 0.5 * N_tri(tri_id, r, s);
    }
}

void RgWedge6Element::dN_tri_dr(int triNodeId, double r, double s, double& dNdr, double& dNds) const
{
    switch (triNodeId) {
        case 0:  // 1-r-s
            dNdr = -1.0; dNds = -1.0;
            break;
        case 1:  // r
            dNdr = 1.0; dNds = 0.0;
            break;
        case 2:  // s
            dNdr = 0.0; dNds = 1.0;
            break;
        default:
            dNdr = 0.0; dNds = 0.0;
    }
}

void RgWedge6Element::evaluateCoordinates(double r, double s, double t,
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

void RgWedge6Element::evaluateJacobian(double r, double s, double t,
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

double RgWedge6Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    // det(J) = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) - ...
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
         - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
         + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

void RgWedge6Element::evaluateJacobianInverse(double r, double s, double t,
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

void RgWedge6Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    // K = integral of B^T * D * B dV using 6-point Gauss quadrature
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
}

void RgWedge6Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgWedge6Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
