#include "RgHex8GeomNLElement.h"
#include "RgFemTypedefs.h"

namespace RgFem {

RgHex8GeomNLElement::RgHex8GeomNLElement() : RgHex8Element() {}

RgHex8GeomNLElement::RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds)
    : RgHex8Element(nodeIds) {}

RgHex8GeomNLElement::RgHex8GeomNLElement(const RgHex8GeomNLElement& other)
    : RgHex8Element(other), m_currentDisplacement(other.m_currentDisplacement) {}

RgHex8GeomNLElement::~RgHex8GeomNLElement() {}

RgElement* RgHex8GeomNLElement::clone() const {
    return new RgHex8GeomNLElement(*this);
}

std::string RgHex8GeomNLElement::typeName() const {
    return "RgHex8GeomNLElement";
}

void RgHex8GeomNLElement::computeStiffness(RgMatrix& Ke, const RgMaterial& material, int integrationOrder) const {
    // Implementation for nonlinear stiffness matrix
    // Typically calls computeTangentStiffnessMatrix internally
    computeTangentStiffnessMatrix(Ke, coords, dN_dx, material, integrationOrder);
}

void RgHex8GeomNLElement::computeTangentStiffness(RgMatrix& Kt, const RgMaterial& material, int integrationOrder) const {
    // Compute the full tangent stiffness matrix
    // Combines material stiffness and geometric (initial stress) stiffness
    // Implementation details depend on material model
}

void RgHex8GeomNLElement::computeStrainStressAtGauss(
    const std::array<std::array<double, 3>, kNodeCount>& coords,
    const std::array<std::array<double, 3>, kNodeCount>& dN_dxi,
    const std::array<double, 3 * kNodeCount>& nodalDisp,
    const RgMaterial& material,
    std::array<double, 6>& strain,
    std::array<double, 6>& stress) const 
{
    // Compute deformation gradient
    std::array<std::array<double, 3>, 3> F;
    computeDeformationGradient(coords, dN_dxi, nodalDisp, F);

    // Compute Green-Lagrange strain
    std::array<std::array<double, 3>, 3> E;
    computeGreenLagrangeStrain(F, E);

    // Flatten 3x3 strain to 6-component vector (Voigt notation)
    strain[0] = E[0][0];          // ε_xx
    strain[1] = E[1][1];          // ε_yy
    strain[2] = E[2][2];          // ε_zz
    strain[3] = 2.0 * E[1][2];    // γ_yz
    strain[4] = 2.0 * E[0][2];    // γ_xz
    strain[5] = 2.0 * E[0][1];    // γ_xy

    // Compute second Piola-Kirchhoff stress
    std::array<std::array<double, 3>, 3> S;
    computeSecondPiolaKirchhoffStress(E, material, S);

    // Flatten 3x3 stress to 6-component vector (Voigt notation)
    stress[0] = S[0][0];
    stress[1] = S[1][1];
    stress[2] = S[2][2];
    stress[3] = S[1][2];
    stress[4] = S[0][2];
    stress[5] = S[0][1];
}

void RgHex8GeomNLElement::computeDeformationGradient(
    const std::array<std::array<double, 3>, kNodeCount>& coords,
    const std::array<std::array<double, 3>, kNodeCount>& dN_dx,
    const std::array<double, 3 * kNodeCount>& nodalDisp,
    std::array<std::array<double, 3>, 3>& F) const 
{
    // Initialize F as identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Compute gradient of displacement field
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double grad_u_ij = 0.0;
            for (int k = 0; k < kNodeCount; ++k) {
                grad_u_ij += nodalDisp[3 * k + j] * dN_dx[k][i];
            }
            F[i][j] += grad_u_ij;
        }
    }
}

void RgHex8GeomNLElement::computeGreenLagrangeStrain(
    const std::array<std::array<double, 3>, 3>& F,
    std::array<std::array<double, 3>, 3>& E) const 
{
    // Compute F^T * F
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

void RgHex8GeomNLElement::computeSecondPiolaKirchhoffStress(
    const std::array<std::array<double, 3>, 3>& E,
    const RgMaterial& material,
    std::array<std::array<double, 3>, 3>& S) const 
{
    // Material-specific implementation (e.g., hyperelasticity)
    // Example: Linearized St. Venant-Kirchhoff material
    double lambda = material.getLameLambda();
    double mu = material.getLameMu();

    double traceE = E[0][0] + E[1][1] + E[2][2];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            S[i][j] = lambda * traceE * ((i == j) ? 1.0 : 0.0) + 2.0 * mu * E[i][j];
        }
    }
}

void RgHex8GeomNLElement::computeTangentStiffnessMatrix(
    RgMatrix& Kt,
    const std::array<std::array<double, 3>, kNodeCount>& coords,
    const std::array<std::array<double, 3>, kNodeCount>& dN_dx,
    const RgMaterial& material,
    int integrationOrder) const 
{
    // Implementation depends on material model and integration scheme
    // Typically combines material stiffness and geometric stiffness
}

} // namespace RgFem