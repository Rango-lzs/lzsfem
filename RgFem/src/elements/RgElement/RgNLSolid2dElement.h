#pragma once

#include "RgSolid2dElement.h"

//! RgNLSolid2dElement - Nonlinear 2D solid element
//! Used for large strain, large displacement analysis (geometrically nonlinear)
//! Features:
//!   - Finite strain formulation
//!   - Updated Lagrangian description
//!   - Deformation gradient and Green-Lagrange strain computation
//!   - Second Piola-Kirchhoff stress
class FEM_EXPORT RgNLSolid2dElement : public RgSolid2dElement
{
public:
    //! default constructor
    RgNLSolid2dElement() = default;

    //! copy constructor
    RgNLSolid2dElement(const RgNLSolid2dElement& el) = default;

    //! assignment operator
    RgNLSolid2dElement& operator=(const RgNLSolid2dElement& el) = default;

    //! destructor
    virtual ~RgNLSolid2dElement() = default;

    //! return element type
    virtual std::string typeName() const override;

    //! Calculate tangent stiffness matrix (nonlinear formulation)
    //! Kt = Km + Kg where Km is material stiffness and Kg is geometric stiffness
    virtual void calculateTangentStiffnessMatrix(RgMatrix& Kt) const override;

    //! Calculate mass matrix (nonlinear, constant in Lagrangian)
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    //! Calculate internal force vector (nonlinear)
    //! F_int = integral of B_nl^T * σ dV in updated Lagrangian formulation
    virtual void calculateInternalForceVector(RgVector& F) const override;

    //! Calculate geometric stiffness (initial stress stiffness)
    //! Kg = integral of B_geo^T * σ * B_geo dV
    virtual void calculateGeometricStiffnessMatrix(RgMatrix& Kg) const override;

    //! Update displacement state for nonlinear iteration
    virtual void updateDisplacementState(const std::vector<double>& displacement) override;

    //! Compute deformation gradient F = I + ∂u/∂X
    virtual void computeDeformationGradient(
        const std::vector<double>& displacement,
        std::array<std::array<double, 2>, 2>& F) const override;

    //! Compute Green-Lagrange strain E = 0.5(C - I) where C = F^T * F
    virtual void computeGreenLagrangeStrain(
        const std::array<std::array<double, 2>, 2>& F,
        std::array<std::array<double, 2>, 2>& E) const override;

    //! Compute Cauchy (true) stress
    virtual void computeCauchyStress(
        const std::array<std::array<double, 2>, 2>& F,
        const RgMaterial& material,
        std::array<std::array<double, 2>, 2>& sigma) const override;
};
