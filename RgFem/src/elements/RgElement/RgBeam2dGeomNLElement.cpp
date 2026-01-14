#include "RgBeam2dGeomNLElement.h"
#include "RgMaterial.h"
#include <cmath>
#include <stdexcept>



RgBeam2dGeomNLElement::RgBeam2dGeomNLElement()
{
}

RgBeam2dGeomNLElement::RgBeam2dGeomNLElement(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgBeam2dGeomNLElement::RgBeam2dGeomNLElement(const RgBeam2dGeomNLElement& other)
    : RgNLBeamElement(other)
{
}

RgBeam2dGeomNLElement& RgBeam2dGeomNLElement::operator=(const RgBeam2dGeomNLElement& other)
{
    if (this != &other) {
        RgNLBeamElement::operator=(other);
    }
    return *this;
}

RgBeam2dGeomNLElement::~RgBeam2dGeomNLElement()
{
}

RgElement* RgBeam2dGeomNLElement::clone() const
{
    return new RgBeam2dGeomNLElement(*this);
}

std::string RgBeam2dGeomNLElement::typeName() const
{
    return "RgBeam2dGeomNLElement";
}

int RgBeam2dGeomNLElement::getNumberOfGaussPoints() const
{
    return 2;  // 2-point Gauss quadrature along beam
}

double RgBeam2dGeomNLElement::evaluateLength() const
{
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    
    double dx = x1[0] - x0[0];
    double dy = x1[1] - x0[1];
    
    return std::sqrt(dx * dx + dy * dy);
}

Vector3d RgBeam2dGeomNLElement::getAxis() const
{
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    
    double length = evaluateLength();
    if (length < 1e-15) {
        return Vector3d(1.0, 0.0, 0.0);
    }
    
    return Vector3d(
        (x1[0] - x0[0]) / length,
        (x1[1] - x0[1]) / length,
        0.0
    );
}

void RgBeam2dGeomNLElement::evaluateShapeFunctions(double r, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    N[0] = N_linear(0, r);  // (1-r)/2
    N[1] = N_linear(1, r);  // (1+r)/2
}

void RgBeam2dGeomNLElement::evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const
{
    dN_dr.resize(kNodeCount);
    dN_dr[0] = dN_linear_dr(0);  // -0.5
    dN_dr[1] = dN_linear_dr(1);  //  0.5
}

double RgBeam2dGeomNLElement::N_linear(int nodeId, double r) const
{
    if (nodeId == 0) {
        return (1.0 - r) / 2.0;
    } else if (nodeId == 1) {
        return (1.0 + r) / 2.0;
    }
    return 0.0;
}

double RgBeam2dGeomNLElement::dN_linear_dr(int nodeId) const
{
    if (nodeId == 0) {
        return -0.5;
    } else if (nodeId == 1) {
        return 0.5;
    }
    return 0.0;
}

void RgBeam2dGeomNLElement::evaluateJacobian(double r, std::array<std::array<double, 2>, 2>& J) const
{
    // For 2D beam in xy-plane
    // J[i][j] = ∂x_i/∂ξ_j where ξ_0 = r
    
    // Initialize
    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
    
    // Only dr column matters for 1D element
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr = dN_linear_dr(node);
        const auto& coord = getNodeCoordinate(node);
        
        J[0][0] += dNdr * coord[0];  // ∂x/∂r
        J[1][0] += dNdr * coord[1];  // ∂y/∂r
    }
}

double RgBeam2dGeomNLElement::evaluateJacobianDeterminant(double r) const
{
    // For 1D element: det(J) = |dx/dr| = magnitude of tangent vector
    std::array<std::array<double, 2>, 2> J;
    evaluateJacobian(r, J);
    
    // For 1D: det = sqrt((dx/dr)^2 + (dy/dr)^2)
    return std::sqrt(J[0][0] * J[0][0] + J[1][0] * J[1][0]);
}

void RgBeam2dGeomNLElement::evaluateJacobianInverse(double r, std::array<std::array<double, 2>, 2>& Jinv) const
{
    std::array<std::array<double, 2>, 2> J;
    evaluateJacobian(r, J);
    
    double det = evaluateJacobianDeterminant(r);
    if (std::abs(det) < 1e-15) {
        Jinv[0][0] = Jinv[1][1] = 1.0;
        Jinv[0][1] = Jinv[1][0] = 0.0;
        return;
    }
    
    // For 1D line element in 2D space:
    // Jinv is computed as pseudo-inverse
    Jinv[0][0] = J[0][0] / (det * det);
    Jinv[0][1] = J[1][0] / (det * det);
    Jinv[1][0] = 0.0;
    Jinv[1][1] = 0.0;
}

void RgBeam2dGeomNLElement::evaluateCoordinates(double r, std::array<double, 3>& coord) const
{
    coord[0] = coord[1] = coord[2] = 0.0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        double N = N_linear(i, r);
        const auto& nodeCoord = getNodeCoordinate(i);
        coord[0] += N * nodeCoord[0];
        coord[1] += N * nodeCoord[1];
        coord[2] += N * nodeCoord[2];
    }
}

void RgBeam2dGeomNLElement::computeLocalAxes(std::array<double, 2>& localX, std::array<double, 2>& localY) const
{
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    
    // Vector along beam
    double dx = x1[0] - x0[0];
    double dy = x1[1] - x0[1];
    double length = std::sqrt(dx * dx + dy * dy);
    
    if (length < 1e-15) {
        // Degenerate element
        localX[0] = 1.0;
        localX[1] = 0.0;
        localY[0] = 0.0;
        localY[1] = 1.0;
        return;
    }
    
    // X-axis along beam
    localX[0] = dx / length;
    localX[1] = dy / length;
    
    // Y-axis perpendicular (rotate X by 90 degrees)
    localY[0] = -localX[1];
    localY[1] = localX[0];
}

void RgBeam2dGeomNLElement::getLocalToGlobalMatrix(std::array<std::array<double, 2>, 2>& R) const
{
    // Rotation matrix from local to global (2D)
    std::array<double, 2> localX, localY;
    computeLocalAxes(localX, localY);
    
    // R[i][j] where i=global, j=local
    R[0][0] = localX[0];  // x_g along x_local
    R[0][1] = localY[0];  // x_g along y_local
    R[1][0] = localX[1];  // y_g along x_local
    R[1][1] = localY[1];  // y_g along y_local
}

void RgBeam2dGeomNLElement::computeDeformationGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                                       std::array<std::array<double, 3>, 3>& F) const
{
    // Initialize F = I
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // For 2D beam: F is 3×3 but we work in 2D plane
    // Get Gauss point coordinate
    double r;
    double weight;
    getGaussPointData(gaussPointIndex, r, weight);
    
    // Compute displacement gradient ∇u
    std::array<std::array<double, 3>, 3> dispGrad;
    computeDisplacementGradient(gaussPointIndex, displacement, dispGrad);
    
    // F = I + ∇u (only in-plane components)
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            F[i][j] += dispGrad[i][j];
        }
    }
}

void RgBeam2dGeomNLElement::computeDisplacementGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                                        std::array<std::array<double, 3>, 3>& dispGrad) const
{
    // Initialize
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dispGrad[i][j] = 0.0;
        }
    }
    
    // Get Gauss point parameter
    double r;
    double weight;
    getGaussPointData(gaussPointIndex, r, weight);
    
    // Get shape function derivatives
    std::vector<double> dN_dr(kNodeCount);
    evaluateShapeFunctions(r, dN_dr);
    evaluateShapeDerivatives(r, dN_dr);
    
    // Get Jacobian inverse
    std::array<std::array<double, 2>, 2> Jinv;
    evaluateJacobianInverse(r, Jinv);
    
    // ∂u/∂x = (∂u/∂r) * (∂r/∂x) = (∂u/∂r) * (∂r/∂x)
    // where x is physical coordinates
    
    // Assemble displacement gradient
    // dispGrad[i][j] = ∂u_i/∂x_j
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr = dN_dr[node];
        
        // Displacement components at this node: [ux, uy, rz]
        double ux = displacement[node * kDofsPerNode + 0];
        double uy = displacement[node * kDofsPerNode + 1];
        
        // Chain rule: ∂u/∂x = (∂u/∂r) * (dr/dx)
        // dr/dx is first column of Jacobian inverse
        dispGrad[0][0] += dNdr * ux * Jinv[0][0];  // ∂ux/∂x
        dispGrad[0][1] += dNdr * ux * Jinv[0][1];  // ∂ux/∂y
        
        dispGrad[1][0] += dNdr * uy * Jinv[0][0];  // ∂uy/∂x
        dispGrad[1][1] += dNdr * uy * Jinv[0][1];  // ∂uy/∂y
    }
}

void RgBeam2dGeomNLElement::getGaussPointData(int gaussPointIndex, double& r, double& weight) const
{
    if (gaussPointIndex == 0) {
        // First Gauss point at -1/sqrt(3)
        r = -1.0 / std::sqrt(3.0);
        weight = 1.0;
    } else if (gaussPointIndex == 1) {
        // Second Gauss point at +1/sqrt(3)
        r = 1.0 / std::sqrt(3.0);
        weight = 1.0;
    } else {
        throw std::out_of_range("Invalid Gauss point index");
    }
}

void RgBeam2dGeomNLElement::calculateStiffnessMatrix(RgMatrix& K) const
{
    int ndofs = kTotalDofs;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // K = K_linear + K_geometric (nonlinear stiffness)
}

void RgBeam2dGeomNLElement::calculateMassMatrix(RgMatrix& M) const
{
    int ndofs = kTotalDofs;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
    // M = integral of N^T * ρ * N dV (consistent mass matrix)
}

void RgBeam2dGeomNLElement::calculateInternalForceVector(RgVector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
    // F_int = integral of B_nl^T * σ dV
}


