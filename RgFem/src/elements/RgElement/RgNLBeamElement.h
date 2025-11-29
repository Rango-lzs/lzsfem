#pragma once

#include "RgBeamElement.h"

namespace RgFem {

/**
 * @class RgNLBeamElement
 * @brief Abstract base class for nonlinear beam elements
 * 
 * Used for geometric nonlinear beam analysis (large deformations, large displacements).
 * Includes deformation gradient computation for hyperelastic and elastoplastic materials.
 * 
 * Theory:
 * - Green-Lagrange strain formulation (material description)
 * - Deformation gradient: F = I + ∇u
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness
 * 
 * Derived classes:
 * - RgBeam3dGeomNLElement: 3D beam with geometric nonlinearity
 */
class FEM_EXPORT RgNLBeamElement : public RgBeamElement
{
public:
    /// Default constructor
    RgNLBeamElement();

    /// Copy constructor
    RgNLBeamElement(const RgNLBeamElement& other);

    /// Assignment operator
    RgNLBeamElement& operator=(const RgNLBeamElement& other);

    /// Destructor
    virtual ~RgNLBeamElement();

    /// Return element type name
    virtual std::string typeName() const override = 0;

    /// Compute deformation gradient F at a Gauss point
    /// @param gaussPointIndex Index of the Gauss point
    /// @param displacement Nodal displacement vector [ux0, uy0, uz0, rx0, ry0, rz0, ux1, uy1, uz1, rx1, ry1, rz1]
    /// @param F Output 3×3 deformation gradient matrix
    virtual void computeDeformationGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                            std::array<std::array<double, 3>, 3>& F) const = 0;

    /// Compute displacement gradient ∇u at a Gauss point
    /// @param gaussPointIndex Index of the Gauss point
    /// @param displacement Nodal displacement vector
    /// @param dispGrad Output 3×3 displacement gradient matrix
    virtual void computeDisplacementGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                             std::array<std::array<double, 3>, 3>& dispGrad) const = 0;

    /// Calculate stiffness matrix (nonlinear formulation with geometric stiffness)
    /// K_L + K_G (linear + geometric stiffness)
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override = 0;

    /// Calculate mass matrix (for dynamic nonlinear analysis)
    /// M = integral of N^T * ρ * N dV
    virtual void calculateMassMatrix(RgMatrix& M) const override = 0;

    /// Calculate internal force vector (nonlinear)
    /// F_int = integral of P^T : F dV (first Piola-Kirchhoff stress)
    virtual void calculateInternalForceVector(RgVector& F) const override = 0;
};

} // namespace RgFem
