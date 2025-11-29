#include "RgNLTri3Element.h"
#include "RgMaterial.h"
#include <cmath>

namespace RgFem {

RgNLTri3Element::RgNLTri3Element()
{
}

RgNLTri3Element::RgNLTri3Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgNLTri3Element::RgNLTri3Element(const RgNLTri3Element& other)
    : RgNLSolid2dElement(other)
{
}

RgNLTri3Element& RgNLTri3Element::operator=(const RgNLTri3Element& other)
{
    if (this != &other) {
        RgNLSolid2dElement::operator=(other);
    }
    return *this;
}

RgNLTri3Element::~RgNLTri3Element()
{
}

RgElement* RgNLTri3Element::clone() const
{
    return new RgNLTri3Element(*this);
}

std::string RgNLTri3Element::typeName() const
{
    return "RgNLTri3Element";
}

double RgNLTri3Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Linear triangular shape functions (barycentric coordinates)
    switch (nodeId) {
        case 0: return L1(r, s);       // 1 - r - s
        case 1: return L2(r);          // r
        case 2: return L3(s);          // s
        default: return 0.0;
    }
}

void RgNLTri3Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                           double& dNdr, double& dNds, double& dNdt) const
{
    dNdt = 0.0;  // No variation in z-direction for 2D element
    
    switch (nodeId) {
        case 0:
            dNdr = -1.0; dNds = -1.0;
            break;
        case 1:
            dNdr = 1.0; dNds = 0.0;
            break;
        case 2:
            dNdr = 0.0; dNds = 1.0;
            break;
        default:
            dNdr = 0.0; dNds = 0.0;
    }
}

void RgNLTri3Element::evaluateCoordinates(double r, double s, double t,
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

void RgNLTri3Element::evaluateJacobian(double r, double s, double t,
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
}

double RgNLTri3Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D triangular element, det(J) = J[0][0]*J[1][1] - J[0][1]*J[1][0]
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

void RgNLTri3Element::evaluateJacobianInverse(double r, double s, double t,
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
    Jinv[2][2] = 1.0;
}

int RgNLTri3Element::getNumberOfGaussPoints() const
{
    return 1;  // 1-point Gauss for constant strain
}

void RgNLTri3Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

void RgNLTri3Element::computeDeformationGradient(int gaussPointIndex,
                                                     const RgVector& displacement,
                                                     std::array<std::array<double, 3>, 3>& F) const
{
    // Initialize identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Compute displacement gradient ∇u
    std::array<std::array<double, 3>, 3> dispGrad;
    computeDisplacementGradient(gaussPointIndex, displacement, dispGrad);
    
    // F = I + ∇u
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] += dispGrad[i][j];
        }
    }
}

void RgNLTri3Element::computeDisplacementGradient(int gaussPointIndex,
                                                      const RgVector& displacement,
                                                      std::array<std::array<double, 3>, 3>& dispGrad) const
{
    // Initialize displacement gradient
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dispGrad[i][j] = 0.0;
        }
    }

    // For 2D triangular element, displacement is constant within element (linear shape functions)
    // ∇u = sum_node (u_node �?∇N_node)
    
    double r = 0.0, s = 0.0;  // Gauss point at centroid for 1-point rule
    
    // Get Jacobian inverse at Gauss point
    std::array<std::array<double, 3>, 3> Jinv;
    evaluateJacobianInverse(r, s, 0.0, Jinv);
    
    for (int node = 0; node < kNodeCount; ++node) {
        // Shape function derivatives in natural coordinates
        double dNdr, dNds, dNdt;
        shapeDerivatives(node, r, s, 0.0, dNdr, dNds, dNdt);
        
        // Physical derivatives: ∂N/∂x = Jinv^T * ∂N/∂�?
        double dNdx = Jinv[0][0] * dNdr + Jinv[1][0] * dNds;
        double dNdy = Jinv[0][1] * dNdr + Jinv[1][1] * dNds;
        
        // Get node displacement (2D: u_x, u_y)
        int dof_x = node * 2;
        int dof_y = node * 2 + 1;
        
        double u_x = (dof_x < displacement.size()) ? displacement[dof_x] : 0.0;
        double u_y = (dof_y < displacement.size()) ? displacement[dof_y] : 0.0;
        
        // Add contribution: ∇u += u �?∇N
        dispGrad[0][0] += u_x * dNdx;  // ∂u_x/∂x
        dispGrad[0][1] += u_x * dNdy;  // ∂u_x/∂y
        dispGrad[1][0] += u_y * dNdx;  // ∂u_y/∂x
        dispGrad[1][1] += u_y * dNdy;  // ∂u_y/∂y
    }
}

void RgNLTri3Element::calculateTangentStiffnessMatrix(RgMatrix& K) const
{
    // K_tan = K_material + K_geometric (initial stress stiffness)
    int ndofs = kNodeCount * 2;  // 2D element
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integrate B^T*D*B over element (material stiffness)
    // Integrate G^T*S*G over element (geometric stiffness)
}

void RgNLTri3Element::calculateGeometricStiffnessMatrix(RgMatrix& Kg) const
{
    // Kg = integral of G^T * S * G dV (initial stress stiffness)
    int ndofs = kNodeCount * 2;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    // Placeholder: actual implementation needed
}

void RgNLTri3Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgNLTri3Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

double RgNLTri3Element::getElementArea() const
{
    return 0.5 * std::abs(evaluateJacobianDeterminant(0.0, 0.0, 0.0));
}

} // namespace RgFem
