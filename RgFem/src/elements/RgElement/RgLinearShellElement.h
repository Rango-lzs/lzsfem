#pragma once

#include "RgElement.h"

namespace RgFem {

/**
 * @class RgLinearShellElement
 * @brief Abstract base class for linear shell elements
 * 
 * Used for small strain, small displacement shell analysis.
 * Geometric nonlinearity is neglected.
 * Combines membrane (in-plane) and bending (out-of-plane) behavior.
 * 
 * Features:
 * - Kirchhoff or Mindlin plate theory
 * - 6 DOF per node (3 translations + 3 rotations)
 * - Thickness property
 * 
 * Derived classes:
 * - RgShell3Element: 3-node triangular shell
 * - RgShell4Element: 4-node bilinear shell
 */
class FEM_EXPORT RgLinearShellElement : public RgElement
{
public:
    /// Default constructor
    RgLinearShellElement();

    /// Copy constructor
    RgLinearShellElement(const RgLinearShellElement& other);

    /// Assignment operator
    RgLinearShellElement& operator=(const RgLinearShellElement& other);

    /// Destructor
    virtual ~RgLinearShellElement();

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

    /// Calculate stiffness matrix (linear formulation)
    /// K = integral of B^T * D * B dV (membrane + bending)
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override = 0;

    /// Calculate mass matrix (linear formulation)
    /// M = integral of N^T * ρ * N dV
    virtual void calculateMassMatrix(RgMatrix& M) const override = 0;

    /// Calculate internal force vector (linear)
    /// F_int = integral of B^T * σ dV
    virtual void calculateInternalForceVector(RgVector& F) const override = 0;
};

} // namespace RgFem
