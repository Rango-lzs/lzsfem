#include "RgNLQuad4Element.h"
#include "RgMaterial.h"
#include <cmath>



RgNLQuad4Element::RgNLQuad4Element()
{
}

RgNLQuad4Element::RgNLQuad4Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgNLQuad4Element::RgNLQuad4Element(const RgNLQuad4Element& other)
    : RgNLSolid2dElement(other)
{
}

RgNLQuad4Element& RgNLQuad4Element::operator=(const RgNLQuad4Element& other)
{
    if (this != &other) {
        RgNLSolid2dElement::operator=(other);
    }
    return *this;
}

RgNLQuad4Element::~RgNLQuad4Element()
{
}

RgElement* RgNLQuad4Element::clone() const
{
    return new RgNLQuad4Element(*this);
}

std::string RgNLQuad4Element::typeName() const
{
    return "RgNLQuad4Element";
}

double RgNLQuad4Element::N(int nodeId, double r, double s) const
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

double RgNLQuad4Element::dNdr(int nodeId, double r, double s) const
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

double RgNLQuad4Element::dNds(int nodeId, double r, double s) const
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

double RgNLQuad4Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    return N(nodeId, r, s);
}

void RgNLQuad4Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                            double& dNdr_out, double& dNds_out, double& dNdt) const
{
    dNdr_out = dNdr(nodeId, r, s);
    dNds_out = dNds(nodeId, r, s);
    dNdt = 0.0;  // No variation in z-direction for 2D element
}

void RgNLQuad4Element::evaluateCoordinates(double r, double s, double t,
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

void RgNLQuad4Element::evaluateJacobian(double r, double s, double t,
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
        
        // Columns: �?∂r, �?∂s
        J[0][0] += dNdr_val * coord[0];
        J[1][0] += dNdr_val * coord[1];
        J[2][0] += dNdr_val * coord[2];
        
        J[0][1] += dNds_val * coord[0];
        J[1][1] += dNds_val * coord[1];
        J[2][1] += dNds_val * coord[2];
    }
}

double RgNLQuad4Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 2D element: det(J) = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

void RgNLQuad4Element::evaluateJacobianInverse(double r, double s, double t,
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

int RgNLQuad4Element::getNumberOfGaussPoints() const
{
    return 4;  // 2×2 Gauss quadrature
}

void RgNLQuad4Element::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

void RgNLQuad4Element::computeDeformationGradient(int gaussPointIndex,
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

void RgNLQuad4Element::computeDisplacementGradient(int gaussPointIndex,
                                                       const RgVector& displacement,
                                                       std::array<std::array<double, 3>, 3>& dispGrad) const
{
    // Initialize displacement gradient
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dispGrad[i][j] = 0.0;
        }
    }

    // Gauss points for 2×2 quadrature
    const double gp[] = { -0.577350269189626, 0.577350269189626 };  // ±1/�?
    
    // Map Gauss point index to r, s coordinates
    int gp_i = gaussPointIndex / 2;
    int gp_j = gaussPointIndex % 2;
    double r = gp[gp_i];
    double s = gp[gp_j];
    
    // Get Jacobian inverse at Gauss point
    std::array<std::array<double, 3>, 3> Jinv;
    evaluateJacobianInverse(r, s, 0.0, Jinv);
    
    for (int node = 0; node < kNodeCount; ++node) {
        // Shape function derivatives in natural coordinates
        double dNdr_val = dNdr(node, r, s);
        double dNds_val = dNds(node, r, s);
        
        // Physical derivatives: ∂N/∂x = Jinv^T * ∂N/∂�?
        double dNdx = Jinv[0][0] * dNdr_val + Jinv[1][0] * dNds_val;
        double dNdy = Jinv[0][1] * dNdr_val + Jinv[1][1] * dNds_val;
        
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

void RgNLQuad4Element::calculateTangentStiffnessMatrix(RgMatrix& K) const
{
    // K_tan = K_material + K_geometric (initial stress stiffness)
    int ndofs = kNodeCount * 2;  // 2D element
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integrate B^T*D*B over element (material stiffness) using 2×2 Gauss quadrature
    // Integrate G^T*S*G over element (geometric stiffness)
}

void RgNLQuad4Element::calculateGeometricStiffnessMatrix(RgMatrix& Kg) const
{
    // Kg = integral of G^T * S * G dV (initial stress stiffness)
    int ndofs = kNodeCount * 2;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    // Placeholder: actual implementation needed
}

void RgNLQuad4Element::calculateMassMatrix(RgMatrix& M) const
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgNLQuad4Element::calculateInternalForceVector(RgVector& F) const
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

double RgNLQuad4Element::getElementArea() const
{
    // Compute area using 2×2 Gauss quadrature
    const double gp[] = { -0.577350269189626, 0.577350269189626 };
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


