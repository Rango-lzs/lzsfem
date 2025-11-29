#include "RgHex8GeomNLElement.h"
#include "RgNLSolid3dElement.h"
#include "RgHex8Element.h"  // For reference to linear shape functions
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>

namespace RgFem {

// ============================================================================
// Gauss Quadrature Data (8-point for hex)
// ============================================================================

const std::array<double, 2> RgHex8GeomNLElement::gaussPoints_1D = {
    -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)
};

const std::array<double, 2> RgHex8GeomNLElement::gaussWeights_1D = {
    1.0, 1.0
};

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8GeomNLElement::RgHex8GeomNLElement()
    : RgNLSolid3dElement()
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds)
    : RgNLSolid3dElement(nodeIds)
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const RgHex8GeomNLElement& other)
    : RgNLSolid3dElement(other)
{
}

RgHex8GeomNLElement::~RgHex8GeomNLElement()
{
}

RgHex8GeomNLElement& RgHex8GeomNLElement::operator=(const RgHex8GeomNLElement& other)
{
    if (this != &other) {
        RgNLSolid3dElement::operator=(other);
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

RgElement* RgHex8GeomNLElement::clone() const {
    return new RgHex8GeomNLElement(*this);
}

std::string RgHex8GeomNLElement::typeName() const {
    return "RgHex8GeomNLElement";
}

// ============================================================================
// Shape Function Methods (8-node linear hex)
// ============================================================================

double RgHex8GeomNLElement::shapeFunction(int nodeId, double r, double s, double t) const
{
    // Linear hexahedral shape functions
    // Node numbering:
    // 0: (-1,-1,-1), 1: (1,-1,-1), 2: (1,1,-1), 3: (-1,1,-1)
    // 4: (-1,-1,1),  5: (1,-1,1),  6: (1,1,1),  7: (-1,1,1)
    
    return N_linear(nodeId, r, s, t);
}

void RgHex8GeomNLElement::shapeDerivatives(int nodeId, double r, double s, double t,
                                            double& dNdr, double& dNds, double& dNdt) const
{
    dN_linear(nodeId, r, s, t, dNdr, dNds, dNdt);
}

double RgHex8GeomNLElement::N_linear(int nodeId, double r, double s, double t) const
{
    double rval, sval, tval;
    
    switch (nodeId) {
        case 0: rval = -1; sval = -1; tval = -1; break;
        case 1: rval =  1; sval = -1; tval = -1; break;
        case 2: rval =  1; sval =  1; tval = -1; break;
        case 3: rval = -1; sval =  1; tval = -1; break;
        case 4: rval = -1; sval = -1; tval =  1; break;
        case 5: rval =  1; sval = -1; tval =  1; break;
        case 6: rval =  1; sval =  1; tval =  1; break;
        case 7: rval = -1; sval =  1; tval =  1; break;
        default: return 0.0;
    }
    
    return 0.125 * (1.0 + r * rval) * (1.0 + s * sval) * (1.0 + t * tval);
}

void RgHex8GeomNLElement::dN_linear(int nodeId, double r, double s, double t,
                                     double& dNdr, double& dNds, double& dNdt) const
{
    double rval, sval, tval;
    
    switch (nodeId) {
        case 0: rval = -1; sval = -1; tval = -1; break;
        case 1: rval =  1; sval = -1; tval = -1; break;
        case 2: rval =  1; sval =  1; tval = -1; break;
        case 3: rval = -1; sval =  1; tval = -1; break;
        case 4: rval = -1; sval = -1; tval =  1; break;
        case 5: rval =  1; sval = -1; tval =  1; break;
        case 6: rval =  1; sval =  1; tval =  1; break;
        case 7: rval = -1; sval =  1; tval =  1; break;
        default: dNdr = dNds = dNdt = 0.0; return;
    }
    
    dNdr = 0.125 * rval * (1.0 + s * sval) * (1.0 + t * tval);
    dNds = 0.125 * (1.0 + r * rval) * sval * (1.0 + t * tval);
    dNdt = 0.125 * (1.0 + r * rval) * (1.0 + s * sval) * tval;
}

// ============================================================================
// Coordinate/Jacobian Methods
// ============================================================================

void RgHex8GeomNLElement::evaluateCoordinates(double r, double s, double t,
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

void RgHex8GeomNLElement::evaluateJacobian(double r, double s, double t,
                                            std::array<std::array<double, 3>, 3>& J) const
{
    // Initialize Jacobian to zero
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            J[i][j] = 0.0;
        }
    }
    
    // J[i][j] = ∂x_i/∂ξ_j where ξ = {r,s,t}
    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr, dNds, dNdt;
        shapeDerivatives(node, r, s, t, dNdr, dNds, dNdt);
        
        const auto& coord = getNodeCoordinate(node);
        
        // ∂x/∂r
        J[0][0] += dNdr * coord[0];
        J[1][0] += dNdr * coord[1];
        J[2][0] += dNdr * coord[2];
        
        // ∂x/∂s
        J[0][1] += dNds * coord[0];
        J[1][1] += dNds * coord[1];
        J[2][1] += dNds * coord[2];
        
        // ∂x/∂t
        J[0][2] += dNdt * coord[0];
        J[1][2] += dNdt * coord[1];
        J[2][2] += dNdt * coord[2];
    }
}

double RgHex8GeomNLElement::evaluateJacobianDeterminant(double r, double s, double t) const
{
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    return matrixDeterminant(J);
}

void RgHex8GeomNLElement::evaluateJacobianInverse(double r, double s, double t,
                                                   std::array<std::array<double, 3>, 3>& Jinv) const
{
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);
    matrixInverse(J, Jinv);
}

// ============================================================================
// Matrix Operations (Determinant and Inverse)
// ============================================================================

double RgHex8GeomNLElement::matrixDeterminant(const std::array<std::array<double, 3>, 3>& A)
{
    // det(A) using cofactor expansion along first row
    double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
               - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
               + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    return det;
}

void RgHex8GeomNLElement::matrixInverse(const std::array<std::array<double, 3>, 3>& A,
                                         std::array<std::array<double, 3>, 3>& Ainv)
{
    double det = matrixDeterminant(A);
    if (std::abs(det) < 1e-15) {
        // Singular matrix, return identity
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Ainv[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return;
    }
    
    // Compute cofactor matrix
    Ainv[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
    Ainv[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det;
    Ainv[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
    
    Ainv[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / det;
    Ainv[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
    Ainv[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det;
    
    Ainv[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
    Ainv[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det;
    Ainv[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;
}

// ============================================================================
// Geometric Nonlinear Methods (from base class)
// ============================================================================

void RgHex8GeomNLElement::computeDeformationGradient(
    int gaussPointIndex,
    const std::vector<double>& displacement,
    std::array<std::array<double, 3>, 3>& F) const
{
    // Initialize F as identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Get natural coordinates of Gauss point
    if (gaussPointIndex >= kGaussPoints) {
        return;  // Invalid Gauss point index
    }
    
    int idx_r = gaussPointIndex / (kGaussPointsPerDir * kGaussPointsPerDir);
    int idx_s = (gaussPointIndex / kGaussPointsPerDir) % kGaussPointsPerDir;
    int idx_t = gaussPointIndex % kGaussPointsPerDir;
    
    double r = gaussPoints_1D[idx_r];
    double s = gaussPoints_1D[idx_s];
    double t = gaussPoints_1D[idx_t];

    // Compute Jacobian inverse at this Gauss point
    std::array<std::array<double, 3>, 3> Jinv;
    evaluateJacobianInverse(r, s, t, Jinv);

    // Compute displacement gradient ∂u/∂X = (∂u/∂x) * J^(-T)
    // First compute physical derivatives of shape functions
    std::array<std::array<double, 3>, 3> J;
    evaluateJacobian(r, s, t, J);

    // Compute displacement gradient in physical space
    std::array<std::array<double, 3>, 3> grad_u;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            grad_u[i][j] = 0.0;
        }
    }

    for (int node = 0; node < kNodeCount; ++node) {
        double dNdr, dNds, dNdt;
        shapeDerivatives(node, r, s, t, dNdr, dNds, dNdt);

        // Transform to physical derivatives: dN/dx = dN/dxi * dxi/dx
        double dN_dx = dNdr * Jinv[0][0] + dNds * Jinv[1][0] + dNdt * Jinv[2][0];
        double dN_dy = dNdr * Jinv[0][1] + dNds * Jinv[1][1] + dNdt * Jinv[2][1];
        double dN_dz = dNdr * Jinv[0][2] + dNds * Jinv[1][2] + dNdt * Jinv[2][2];

        // Add contribution to displacement gradient
        if (3 * node < displacement.size()) {
            double u_x = displacement[3 * node];
            double u_y = (3 * node + 1 < displacement.size()) ? displacement[3 * node + 1] : 0.0;
            double u_z = (3 * node + 2 < displacement.size()) ? displacement[3 * node + 2] : 0.0;

            grad_u[0][0] += u_x * dN_dx;
            grad_u[0][1] += u_x * dN_dy;
            grad_u[0][2] += u_x * dN_dz;

            grad_u[1][0] += u_y * dN_dx;
            grad_u[1][1] += u_y * dN_dy;
            grad_u[1][2] += u_y * dN_dz;

            grad_u[2][0] += u_z * dN_dx;
            grad_u[2][1] += u_z * dN_dy;
            grad_u[2][2] += u_z * dN_dz;
        }
    }

    // F = I + ∇u
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = ((i == j) ? 1.0 : 0.0) + grad_u[i][j];
        }
    }
}

void RgHex8GeomNLElement::computeGreenLagrangeStrain(
    const std::array<std::array<double, 3>, 3>& F,
    std::array<std::array<double, 3>, 3>& E) const
{
    // Compute F^T * F (Right Cauchy-Green tensor)
    std::array<std::array<double, 3>, 3> FTF;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            FTF[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                FTF[i][j] += F[k][i] * F[k][j];
            }
        }
    }

    // Compute E = 0.5*(F^T*F - I)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            E[i][j] = 0.5 * (FTF[i][j] - ((i == j) ? 1.0 : 0.0));
        }
    }
}

void RgHex8GeomNLElement::computeCauchyStress(
    const std::array<std::array<double, 3>, 3>& F,
    const RgMaterial& material,
    std::array<std::array<double, 3>, 3>& sigma) const
{
    // Compute Green-Lagrange strain
    std::array<std::array<double, 3>, 3> E;
    computeGreenLagrangeStrain(F, E);

    // For now, use simplified hyperelastic material (St. Venant-Kirchhoff)
    // σ = (1/det(F)) * F * S * F^T where S = λ*tr(E)*I + 2*μ*E
    
    double lambda = material.getLameLambda();
    double mu = material.getLameMu();
    
    // Compute second Piola-Kirchhoff stress
    double traceE = E[0][0] + E[1][1] + E[2][2];
    
    std::array<std::array<double, 3>, 3> S;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            S[i][j] = 2.0 * mu * E[i][j];
            if (i == j) {
                S[i][j] += lambda * traceE;
            }
        }
    }

    // Compute Cauchy stress σ = (1/det(F)) * F * S * F^T
    double detF = F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1])
                - F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0])
                + F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]);

    if (std::abs(detF) < 1e-15) {
        // Singular deformation, return zero stress
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                sigma[i][j] = 0.0;
            }
        }
        return;
    }

    // F * S
    std::array<std::array<double, 3>, 3> FS;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            FS[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                FS[i][j] += F[i][k] * S[k][j];
            }
        }
    }

    // (F*S)*F^T / det(F)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            sigma[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                sigma[i][j] += FS[i][k] * F[j][k];
            }
            sigma[i][j] /= detF;
        }
    }
}

// ============================================================================
// Matrix Assembly Methods (Stub implementations)
// ============================================================================

void RgHex8GeomNLElement::calculateTangentStiffnessMatrix(RgMatrix& Kt) const
{
    // Tangent stiffness for Newton-Raphson: Kt = Km + Kg
    // Km: material (constitutive) stiffness
    // Kg: geometric (initial stress) stiffness
    int ndofs = kNodeCount * 3;
    Kt.resize(ndofs, ndofs);
    Kt.zero();
    
    // Placeholder: would compute combined stiffness
}

void RgHex8GeomNLElement::calculateMassMatrix(RgMatrix& M) const
{
    // Mass matrix (constant in Lagrangian description)
    // M = integral of N^T * ρ * N dV
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: would compute mass matrix
}

void RgHex8GeomNLElement::calculateInternalForceVector(RgVector& F) const
{
    // Internal force: F_int = integral of B^T * σ dV
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    // Placeholder: would compute stress-driven internal forces
}

void RgHex8GeomNLElement::calculateGeometricStiffnessMatrix(RgMatrix& Kg) const
{
    // Geometric stiffness: Kg from initial stress effects
    int ndofs = kNodeCount * 3;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    // Placeholder: would include stress-stiffness coupling
}

// ============================================================================
// State Management
// ============================================================================

void RgHex8GeomNLElement::updateDisplacementState(const std::vector<double>& displacement)
{
    // Store updated displacement for next iteration
    m_currentDisplacement = displacement;
}

} // namespace RgFem