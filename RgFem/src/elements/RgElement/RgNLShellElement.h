#pragma once

#include "RgShellElement.h"
#include <array>
#include <vector>



/**
 * @class RgNLShellElement
 * @brief Abstract base class for nonlinear shell elements
 * 
 * Used for geometric nonlinear shell analysis (large deformations, large displacements).
 * Combines membrane and bending behavior with geometric stiffness effects.
 * Includes deformation gradient computation for hyperelastic and elastoplastic materials.
 * 
 * Theory:
 * - Green-Lagrange strain formulation (material description)
 * - Deformation gradient: F = I + ∇u
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness
 * - 6 DOF per node (3 translations + 3 rotations)
 * 
 * Derived classes:
 * - RgShell3GeomNLElement: 3-node triangular shell with geometric nonlinearity
 * - RgShell4GeomNLElement: 4-node bilinear shell with geometric nonlinearity
 */
class FEM_EXPORT RgNLShellElement : public RgShellElement
{
public:
    /// Default constructor
    RgNLShellElement();

    /// Copy constructor
    RgNLShellElement(const RgNLShellElement& other);

    /// Assignment operator
    RgNLShellElement& operator=(const RgNLShellElement& other);

    /// Destructor
    virtual ~RgNLShellElement();

    /// Return element type name
    virtual std::string typeName() const override = 0;

    /// Get shell thickness
    virtual double getShellThickness() const = 0;

    /// Set shell thickness
    virtual void setShellThickness(double thickness) = 0;

    /// Get element area
    virtual double getElementArea() const = 0;

    /// Get normal vector at parametric coordinates (r, s)
    virtual void getShellNormal(double r, double s, std::array<double, 3>& normal) const = 0;

    /// Compute deformation gradient F at a Gauss point
    /// @param gaussPointIndex Index of the Gauss point
    /// @param displacement Nodal displacement vector
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


