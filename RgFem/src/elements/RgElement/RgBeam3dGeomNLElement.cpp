#include "RgBeam3dGeomNLElement.h"
#include "RgMaterial.h"
#include <cmath>

namespace RgFem {

RgBeam3dGeomNLElement::RgBeam3dGeomNLElement()
{
}

RgBeam3dGeomNLElement::RgBeam3dGeomNLElement(const std::array<int, kNodeCount>& nodeIds)
{
    setNodeIds(nodeIds);
}

RgBeam3dGeomNLElement::RgBeam3dGeomNLElement(const RgBeam3dGeomNLElement& other)
    : RgBeamElement(other)
{
}

RgBeam3dGeomNLElement& RgBeam3dGeomNLElement::operator=(const RgBeam3dGeomNLElement& other)
{
    if (this != &other) {
        RgBeamElement::operator=(other);
    }
    return *this;
}

RgBeam3dGeomNLElement::~RgBeam3dGeomNLElement()
{
}

RgElement* RgBeam3dGeomNLElement::clone() const
{
    return new RgBeam3dGeomNLElement(*this);
}

std::string RgBeam3dGeomNLElement::typeName() const
{
    return "RgBeam3dGeomNLElement";
}

double RgBeam3dGeomNLElement::N_cubic(int nodeId, double r) const
{
    // Cubic Hermite shape functions for 2-node beam element
    // r in [-1, 1]: r = -1 at node 0, r = +1 at node 1
    // Transformed to local coordinate: xi = (r+1)/2, xi in [0,1]
    
    double xi = (r + 1.0) / 2.0;
    
    if (nodeId == 0) {
        // Node 0 shape function: N0 = 1 - 3*xi^2 + 2*xi^3
        return 1.0 - 3.0 * xi * xi + 2.0 * xi * xi * xi;
    } else if (nodeId == 1) {
        // Node 1 shape function: N1 = 3*xi^2 - 2*xi^3
        return 3.0 * xi * xi - 2.0 * xi * xi * xi;
    }
    return 0.0;
}

double RgBeam3dGeomNLElement::dN_cubic_dr(int nodeId, double r) const
{
    // Derivative with respect to r: d/dr = (2/L) * d/dxi
    // For normalized [-1, 1] domain: d/dr = 2 * d/dxi
    
    double xi = (r + 1.0) / 2.0;
    
    if (nodeId == 0) {
        // dN0/dxi = -6*xi + 6*xi^2
        double dNdxi = -6.0 * xi + 6.0 * xi * xi;
        return 2.0 * dNdxi;
    } else if (nodeId == 1) {
        // dN1/dxi = 6*xi - 6*xi^2
        double dNdxi = 6.0 * xi - 6.0 * xi * xi;
        return 2.0 * dNdxi;
    }
    return 0.0;
}

double RgBeam3dGeomNLElement::shapeFunction(int nodeId, double r, double s, double t) const
{
    // For beam element, only variation along r-axis
    return N_cubic(nodeId, r);
}

void RgBeam3dGeomNLElement::shapeDerivatives(int nodeId, double r, double s, double t,
                                            double& dNdr, double& dNds, double& dNdt) const
{
    dNdr = dN_cubic_dr(nodeId, r);
    dNds = 0.0;  // No variation in s or t
    dNdt = 0.0;
}

void RgBeam3dGeomNLElement::evaluateCoordinates(double r, double s, double t,
                                               std::array<double, 3>& coord) const
{
    coord[0] = coord[1] = coord[2] = 0.0;
    
    for (int i = 0; i < kNodeCount; ++i) {
        double N = N_cubic(i, r);
        const auto& nodeCoord = getNodeCoordinate(i);
        coord[0] += N * nodeCoord[0];
        coord[1] += N * nodeCoord[1];
        coord[2] += N * nodeCoord[2];
    }
}

void RgBeam3dGeomNLElement::evaluateJacobian(double r, double s, double t,
                                            std::array<std::array<double, 3>, 3>& J) const
{
    // Initialize Jacobian
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            J[i][j] = 0.0;
        }
    }

    // For beam, only variation along r-axis
    // J = dx/dr
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr = dN_cubic_dr(node, r);
        const auto& coord = getNodeCoordinate(node);
        
        J[0][0] += dNdr * coord[0];
        J[1][0] += dNdr * coord[1];
        J[2][0] += dNdr * coord[2];
    }
}

double RgBeam3dGeomNLElement::evaluateJacobianDeterminant(double r, double s, double t) const
{
    // For 1D beam element, det(J) = |dx/dr| (length of element)
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    double dx_dr = std::sqrt(J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]);
    return dx_dr;
}

void RgBeam3dGeomNLElement::evaluateJacobianInverse(double r, double s, double t,
                                                   std::array<std::array<double, 3>, 3>& Jinv) const
{
    // For 1D beam element, inverse is simpler
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Jinv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    
    double det = evaluateJacobianDeterminant(r, s, t);
    if (std::abs(det) < 1e-15) {
        return;
    }
    
    // Jinv[i][0] = J[0][i] / det (inverse for 1D)
    Jinv[0][0] = J[0][0] / (det * det);
    Jinv[1][0] = J[1][0] / (det * det);
    Jinv[2][0] = J[2][0] / (det * det);
}

int RgBeam3dGeomNLElement::getNumberOfGaussPoints() const
{
    return 2;  // 2-point Gauss quadrature along beam axis
}

void RgBeam3dGeomNLElement::initTraits()
{
    // Initialize element traits (Gauss points, shape function derivatives)
    // To be implemented
}

double RgBeam3dGeomNLElement::getBeamLength() const
{
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    
    double dx = x1[0] - x0[0];
    double dy = x1[1] - x0[1];
    double dz = x1[2] - x0[2];
    
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void RgBeam3dGeomNLElement::computeLocalAxes(std::array<double, 3>& xAxis,
                                            std::array<double, 3>& yAxis,
                                            std::array<double, 3>& zAxis) const
{
    const auto& x0 = getNodeCoordinate(0);
    const auto& x1 = getNodeCoordinate(1);
    
    // x-axis along beam
    double L = getBeamLength();
    if (L > 1e-15) {
        xAxis[0] = (x1[0] - x0[0]) / L;
        xAxis[1] = (x1[1] - x0[1]) / L;
        xAxis[2] = (x1[2] - x0[2]) / L;
    } else {
        xAxis[0] = 1.0; xAxis[1] = 0.0; xAxis[2] = 0.0;
    }
    
    // y-axis perpendicular to x (choose z-direction if x is horizontal)
    if (std::abs(xAxis[2]) < 0.999) {
        yAxis[0] = -xAxis[0] * xAxis[1] / (1.0 - xAxis[2]);
        yAxis[1] = 1.0 - xAxis[1] * xAxis[1] / (1.0 - xAxis[2]);
        yAxis[2] = 0.0;
    } else {
        yAxis[0] = 1.0;
        yAxis[1] = 0.0;
        yAxis[2] = 0.0;
    }
    
    // Normalize y-axis
    double ySqrt = std::sqrt(yAxis[0] * yAxis[0] + yAxis[1] * yAxis[1] + yAxis[2] * yAxis[2]);
    if (ySqrt > 1e-15) {
        yAxis[0] /= ySqrt;
        yAxis[1] /= ySqrt;
        yAxis[2] /= ySqrt;
    }
    
    // z-axis = x × y (cross product)
    zAxis[0] = xAxis[1] * yAxis[2] - xAxis[2] * yAxis[1];
    zAxis[1] = xAxis[2] * yAxis[0] - xAxis[0] * yAxis[2];
    zAxis[2] = xAxis[0] * yAxis[1] - xAxis[1] * yAxis[0];
}

void RgBeam3dGeomNLElement::getLocalXAxis(std::array<double, 3>& xAxis) const
{
    std::array<double, 3> yAxis, zAxis;
    computeLocalAxes(xAxis, yAxis, zAxis);
}

void RgBeam3dGeomNLElement::getLocalYAxis(std::array<double, 3>& yAxis) const
{
    std::array<double, 3> xAxis, zAxis;
    computeLocalAxes(xAxis, yAxis, zAxis);
}

void RgBeam3dGeomNLElement::getLocalZAxis(std::array<double, 3>& zAxis) const
{
    std::array<double, 3> xAxis, yAxis;
    computeLocalAxes(xAxis, yAxis, zAxis);
}

void RgBeam3dGeomNLElement::getLocalToGlobalMatrix(std::array<std::array<double, 3>, 3>& T) const
{
    std::array<double, 3> xAxis, yAxis, zAxis;
    computeLocalAxes(xAxis, yAxis, zAxis);
    
    // Transformation matrix: rows are local axes
    for (int i = 0; i < 3; ++i) {
        T[0][i] = xAxis[i];
        T[1][i] = yAxis[i];
        T[2][i] = zAxis[i];
    }
}

void RgBeam3dGeomNLElement::computeDeformationGradient(int gaussPointIndex,
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

void RgBeam3dGeomNLElement::computeDisplacementGradient(int gaussPointIndex,
                                                       const RgVector& displacement,
                                                       std::array<std::array<double, 3>, 3>& dispGrad) const
{
    // Initialize displacement gradient
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dispGrad[i][j] = 0.0;
        }
    }

    // Gauss points for 2-point quadrature: ±1/√3
    const double gp[] = { -0.577350269189626, 0.577350269189626 };
    double r = gp[gaussPointIndex];
    
    // ∇u = du/dr for beam (1D analysis)
    // du/dr = sum_node (du_node/dr)
    
    double dudx = 0.0, dvdx = 0.0, dwdx = 0.0;
    
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr = dN_cubic_dr(node, r);
        
        // Get node displacement (first 3 DOF per node: u, v, w)
        int dof_u = node * kDofsPerNode;
        int dof_v = node * kDofsPerNode + 1;
        int dof_w = node * kDofsPerNode + 2;
        
        double u = (dof_u < displacement.size()) ? displacement[dof_u] : 0.0;
        double v = (dof_v < displacement.size()) ? displacement[dof_v] : 0.0;
        double w = (dof_w < displacement.size()) ? displacement[dof_w] : 0.0;
        
        dudx += u * dNdr;
        dvdx += v * dNdr;
        dwdx += w * dNdr;
    }
    
    // Account for Jacobian: du/dx = du/dr * dr/dx
    double det = evaluateJacobianDeterminant(r, 0.0, 0.0);
    if (std::abs(det) > 1e-15) {
        dudx /= det;
        dvdx /= det;
        dwdx /= det;
    }
    
    // Fill displacement gradient
    dispGrad[0][0] = dudx;
    dispGrad[1][1] = dvdx;
    dispGrad[2][2] = dwdx;
}

void RgBeam3dGeomNLElement::calculateTangentStiffnessMatrix(RgMatrix& K) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
}

void RgBeam3dGeomNLElement::calculateGeometricStiffnessMatrix(RgMatrix& Kg) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    // Placeholder: actual implementation needed
}

void RgBeam3dGeomNLElement::calculateMassMatrix(RgMatrix& M) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgBeam3dGeomNLElement::calculateInternalForceVector(RgVector& F) const
{
    int ndofs = kNodeCount * kDofsPerNode;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: actual implementation needed
}

} // namespace RgFem
