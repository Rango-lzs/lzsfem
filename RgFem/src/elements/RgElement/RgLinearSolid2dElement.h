#pragma once

#include "RgSolid2dElement.h"

//! RgLinearSolid2dElement - Linear 2D solid element
//! Used for small strain, small displacement analysis (geometrically linear)
class FEM_EXPORT RgLinearSolid2dElement : public RgSolid2dElement
{
public:
    //! default constructor
    RgLinearSolid2dElement() = default;

    //! copy constructor
    RgLinearSolid2dElement(const RgLinearSolid2dElement& el) = default;

    //! assignment operator
    RgLinearSolid2dElement& operator=(const RgLinearSolid2dElement& el) = default;

    //! destructor
    virtual ~RgLinearSolid2dElement() = default;

    //! return element type
    virtual std::string typeName() const override;

    //! Calculate stiffness matrix (linear formulation)
    //! K = integral of B^T * D * B dV
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;

    //! Calculate mass matrix (linear formulation)
    //! M = integral of N^T * ρ * N dV
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    //! Calculate internal force vector (linear)
    //! F_int = integral of B^T * σ dV
    virtual void calculateInternalForceVector(RgVector& F) const override;
};
