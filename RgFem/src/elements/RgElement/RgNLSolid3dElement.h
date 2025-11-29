#pragma once

#include "RgSolid3dElement.h"

//! RgNLSolid3dElement - Nonlinear 3D solid element
//! Used for large strain, large displacement analysis (geometrically nonlinear)
//! Features:
//!   - Finite strain formulation
//!   - Updated Lagrangian description
//!   - Deformation gradient F = I + ∂u/∂X
//!   - Green-Lagrange strain E = 0.5(C - I)
//!   - Cauchy (true) stress computation
//!   - Geometric stiffness from initial stress effects
//!   - 8-point Gauss quadrature for 3D elements
class FEM_EXPORT RgNLSolid3dElement : public RgSolid3dElement
{
public:
    //! default constructor
    RgNLSolid3dElement() = default;

    //! copy constructor
    RgNLSolid3dElement(const RgNLSolid3dElement& el) = default;

    //! assignment operator
    RgNLSolid3dElement& operator=(const RgNLSolid3dElement& el) = default;

    //! destructor
    virtual ~RgNLSolid3dElement() = default;

    //! return element type
    virtual std::string typeName() const override;

    //! Calculate tangent stiffness matrix (nonlinear formulation)
    //! Kt = Km + Kg where:
    //!   - Km = material stiffness (constitutive tangent)
    //!   - Kg = geometric stiffness (initial stress effect)
    virtual void calculateTangentStiffnessMatrix(RgMatrix& Kt) const override;

    //! Calculate mass matrix (nonlinear, constant in Lagrangian)
    //! M = integral of N^T * ρ * N dV (reference configuration)
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    //! Calculate internal force vector (nonlinear)
    //! F_int = integral of B_nl^T * σ dV in updated Lagrangian formulation
    virtual void calculateInternalForceVector(RgVector& F) const override;

    //! Calculate geometric stiffness (initial stress stiffness)
    //! Kg = integral of B_geo^T * σ * B_geo dV
    virtual void calculateGeometricStiffnessMatrix(RgMatrix& Kg) const override;

    //! Update displacement state for nonlinear iteration
    //! Caches deformation gradients and stresses at all Gauss points
    virtual void updateDisplacementState(const std::vector<double>& displacement) override;

    //! Compute deformation gradient F = I + ∂u/∂X
    //! Input: displacement vector
    //! Output: 3×3 deformation gradient tensor at specified Gauss point
    virtual void computeDeformationGradient(
        int gaussPointIndex,
        const std::vector<double>& displacement,
        std::array<std::array<double, 3>, 3>& F) const override;

    //! Compute Green-Lagrange strain E = 0.5(C - I) where C = F^T * F
    virtual void computeGreenLagrangeStrain(
        const std::array<std::array<double, 3>, 3>& F,
        std::array<std::array<double, 3>, 3>& E) const override;

    //! Compute Cauchy (true) stress from deformation and material
    //! σ = (1/det(F)) * F * S * F^T where S is Second Piola-Kirchhoff stress
    virtual void computeCauchyStress(
        const std::array<std::array<double, 3>, 3>& F,
        const RgMaterial& material,
        std::array<std::array<double, 3>, 3>& sigma) const override;

    //! Compute material stiffness matrix (tangent constitutive matrix)
    virtual void computeMaterialStiffnessMatrix(
        int gaussPointIndex,
        const RgMaterial& material,
        RgMatrix& Cm) const override;

    //! Get cached deformation gradient at Gauss point
    virtual const std::array<std::array<double, 3>, 3>& getDeformationGradient(int gaussPointIndex) const override;

    //! Get cached Cauchy stress at Gauss point
    virtual const std::array<std::array<double, 3>, 3>& getCauchyStress(int gaussPointIndex) const override;
};
