#include "RgShell4Element.h"
#include "RgMaterial.h"
#include <cmath>

namespace RgFem {

RgShell4Element::RgShell4Element()
{
}

RgShell4Element::RgShell4Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgShell4Element::RgShell4Element(const RgShell4Element& other)
    : RgElement(other), m_thickness(other.m_thickness)
{
}

RgShell4Element& RgShell4Element::operator=(const RgShell4Element& other)
{
    if (this != &other) {
        RgElement::operator=(other);
        m_thickness = other.m_thickness;
    }
    return *this;
}

RgShell4Element::~RgShell4Element()
{
}

RgElement* RgShell4Element::clone() const
{
    return new RgShell4Element(*this);
}

std::string RgShell4Element::typeName() const
{
    return "RgShell4Element";
}

double RgShell4Element::N(int nodeId, double r, double s) const
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

double RgShell4Element::dNdr(int nodeId, double r, double s) const
{
    switch (nodeId) {
        case 0: return -0.25 * (1.0 - s);
        case 1: return  0.25 * (1.0 - s);
        case 2: return  0.25 * (1.0 + s);
        case 3: return -0.25 * (1.0 + s);
        default: return 0.0;
    }
}

double RgShell4Element::dNds(int nodeId, double r, double s) const
{
    switch (nodeId) {
        case 0: return -0.25 * (1.0 - r);
        case 1: return -0.25 * (1.0 + r);
        case 2: return  0.25 * (1.0 + r);
        case 3: return  0.25 * (1.0 - r);
        default: return 0.0;
    }
}

double RgShell4Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    return N(nodeId, r, s);
}

void RgShell4Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                       double& dNdr_out, double& dNds_out, double& dNdt) const
{
    dNdr_out = dNdr(nodeId, r, s);
    dNds_out = dNds(nodeId, r, s);
    dNdt = 0.0;
}

void RgShell4Element::evaluateCoordinates(double r, double s, double t,
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

void RgShell4Element::evaluateJacobian(double r, double s, double t,
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
        
        J[0][0] += dNdr_val * coord[0];
        J[1][0] += dNdr_val * coord[1];
        J[2][0] += dNdr_val * coord[2];
        
        J[0][1] += dNds_val * coord[0];
        J[1][1] += dNds_val * coord[1];
        J[2][1] += dNds_val * coord[2];
    }
}

double RgShell4Element::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For shell element (2D surface): det(J) = |∂x/∂r × ∂x/∂s|
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    // Cross product of columns 0 and 1
    double cx = J[1][0] * J[2][1] - J[2][0] * J[1][1];
    double cy = J[2][0] * J[0][1] - J[0][0] * J[2][1];
    double cz = J[0][0] * J[1][1] - J[1][0] * J[0][1];
    
    return std::sqrt(cx * cx + cy * cy + cz * cz);
}

void RgShell4Element::evaluateJacobianInverse(double r, double s, double t,
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
    
    // For 2D shell in 3D space: compute pseudo-inverse
    // (J^T J)^-1 J^T
    double det2 = det * det;
    
    // Simplified 2×2 inverse of the Jacobian metric
    Jinv[0][0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) / det2;
    Jinv[0][1] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) / det2;
    Jinv[1][0] = (J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0]) / det2;
    Jinv[1][1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) / det2;
}

int RgShell4Element::getNumberOfGaussPoints() const
{
    return 4;  // 2×2 in-plane quadrature
}

void RgShell4Element::initTraits()
{
    // Initialize element traits
    // To be implemented
}

void RgShell4Element::getShellNormal(double r, double s, std::array<double, 3>& normal) const
{
    // Normal = ∂x/∂r × ∂x/∂s
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, 0.0, J);
    
    // Cross product of columns 0 and 1
    normal[0] = J[1][0] * J[2][1] - J[2][0] * J[1][1];
    normal[1] = J[2][0] * J[0][1] - J[0][0] * J[2][1];
    normal[2] = J[0][0] * J[1][1] - J[1][0] * J[0][1];
    
    // Normalize
    double mag = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (mag > 1e-15) {
        normal[0] /= mag;
        normal[1] /= mag;
        normal[2] /= mag;
    }
}

double RgShell4Element::getShellThickness() const
{
    return m_thickness;
}

void RgShell4Element::setShellThickness(double thickness)
{
    m_thickness = thickness;
}

double RgShell4Element::getElementArea() const
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
    
    return area;
}

void RgShell4Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integrate membrane + bending stiffness over element and thickness
}

void RgShell4Element::calculateMassMatrix(RgMatrix& M) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgShell4Element::calculateInternalForceVector(RgVector& F) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
