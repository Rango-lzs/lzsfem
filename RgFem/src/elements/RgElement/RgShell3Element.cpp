#include "RgShell3Element.h"
#include "RgMaterial.h"
#include <cmath>

namespace RgFem {

RgShell3Element::RgShell3Element()
{
}

RgShell3Element::RgShell3Element(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgShell3Element::RgShell3Element(const RgShell3Element& other)
    : RgElement(other), m_thickness(other.m_thickness)
{
}

RgShell3Element& RgShell3Element::operator=(const RgShell3Element& other)
{
    if (this != &other) {
        RgElement::operator=(other);
        m_thickness = other.m_thickness;
    }
    return *this;
}

RgShell3Element::~RgShell3Element()
{
}

RgElement* RgShell3Element::clone() const
{
    return new RgShell3Element(*this);
}

std::string RgShell3Element::typeName() const
{
    return "RgShell3Element";
}

double RgShell3Element::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Linear triangular shape functions (barycentric coordinates)
    switch (nodeId) {
        case 0: return L1(r, s);       // 1 - r - s
        case 1: return L2(r);          // r
        case 2: return L3(s);          // s
        default: return 0.0;
    }
}

void RgShell3Element::shapeDerivatives(int nodeId, double r, double s, double t,
                                       double& dNdr, double& dNds, double& dNdt) const
{
    dNdt = 0.0;
    
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

void RgShell3Element::evaluateCoordinates(double r, double s, double t,
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

void RgShell3Element::evaluateJacobian(double r, double s, double t,
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

double RgShell3Element::evaluateJacobianDeterminant(double r, double s, double t) const
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

void RgShell3Element::evaluateJacobianInverse(double r, double s, double t,
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
    double det2 = det * det;
    
    Jinv[0][0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]) / det2;
    Jinv[0][1] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]) / det2;
    Jinv[1][0] = (J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0]) / det2;
    Jinv[1][1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]) / det2;
}

int RgShell3Element::getNumberOfGaussPoints() const
{
    return 1;  // 1-point in-plane quadrature (constant strain triangle)
}

void RgShell3Element::initTraits()
{
    // Initialize element traits
    // To be implemented
}

void RgShell3Element::getShellNormal(double r, double s, std::array<double, 3>& normal) const
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

double RgShell3Element::getShellThickness() const
{
    return m_thickness;
}

void RgShell3Element::setShellThickness(double thickness)
{
    m_thickness = thickness;
}

double RgShell3Element::getElementArea() const
{
    // For constant strain triangle, area is constant
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    const auto& x2 = getNodeCoordinate(2);
    
    // Vectors from node 0 to nodes 1 and 2
    double v1_x = x1[0] - x0[0];
    double v1_y = x1[1] - x0[1];
    double v1_z = x1[2] - x0[2];
    double v2_x = x2[0] - x0[0];
    double v2_y = x2[1] - x0[1];
    double v2_z = x2[2] - x0[2];
    
    // Cross product
    double cx = v1_y * v2_z - v1_z * v2_y;
    double cy = v1_z * v2_x - v1_x * v2_z;
    double cz = v1_x * v2_y - v1_y * v2_x;
    
    double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
    return area;
}

void RgShell3Element::calculateStiffnessMatrix(RgMatrix& K) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integrate membrane + bending stiffness over element and thickness
}

void RgShell3Element::calculateMassMatrix(RgMatrix& M) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgShell3Element::calculateInternalForceVector(RgVector& F) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
