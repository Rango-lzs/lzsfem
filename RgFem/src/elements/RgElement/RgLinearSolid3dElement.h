#pragma once

#include "RgSolid3dElement.h"

//! RgLinearSolid3dElement - Linear 3D solid element
//! Used for small strain, small displacement analysis (geometrically linear)
//! Features:
//!   - Linear strain-displacement relationship
//!   - Standard B-matrix formulation
//!   - Material stiffness integration via Gauss quadrature
//!   - Compatible with isotropic and anisotropic materials
class FEM_EXPORT RgLinearSolid3dElement : public RgSolid3dElement
{
public:
    //! default constructor
    RgLinearSolid3dElement() = default;

    //! copy constructor
    RgLinearSolid3dElement(const RgLinearSolid3dElement& el) = default;

    //! assignment operator
    RgLinearSolid3dElement& operator=(const RgLinearSolid3dElement& el) = default;

    //! destructor
    virtual ~RgLinearSolid3dElement() = default;

    //! return element type
    virtual std::string typeName() const override;

    //! Calculate stiffness matrix (linear formulation)
    //! K = integral of B^T * D * B dV
    //! Using Gauss quadrature integration
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;

    //! Calculate mass matrix (linear formulation)
    //! M = integral of N^T * ρ * N dV
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    //! Calculate internal force vector (linear)
    //! F_int = integral of B^T * σ dV
    virtual void calculateInternalForceVector(RgVector& F) const override;

    //! Compute strain-displacement matrix B at a Gauss point
    //! Relates nodal displacements to strains: ε = B * u
    virtual void computeStrainDisplacementMatrix(
        const std::vector<double>& dN_dx,
        const std::vector<double>& dN_dy,
        const std::vector<double>& dN_dz,
        RgMatrix& B) const override;

    //! Compute physical derivatives of shape functions at natural coordinates
    virtual void computePhysicalDerivatives(
        double r, double s, double t,
        std::vector<double>& dN_dx,
        std::vector<double>& dN_dy,
        std::vector<double>& dN_dz,
        double& jacobianDet) const override;
};
