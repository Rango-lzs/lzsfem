#pragma once

#include "RgBeamElement.h"

namespace RgFem {

/**
 * @class RgLinearBeamElement
 * @brief Abstract base class for linear beam elements
 * 
 * Used for small strain, small displacement, small rotation beam analysis.
 * Geometric nonlinearity is neglected.
 * 
 * Derived classes:
 * - RgBeam2dElement: 2D planar beam
 * - RgBeam3dElement: 3D spatial beam
 */
class FEM_EXPORT RgLinearBeamElement : public RgBeamElement
{
public:
    /// Default constructor
    RgLinearBeamElement();

    /// Copy constructor
    RgLinearBeamElement(const RgLinearBeamElement& other);

    /// Assignment operator
    RgLinearBeamElement& operator=(const RgLinearBeamElement& other);

    /// Destructor
    virtual ~RgLinearBeamElement();

    /// Return element type name
    virtual std::string typeName() const override = 0;

    /// Calculate stiffness matrix (linear formulation)
    /// K = integral of B^T * D * B dV
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override = 0;

    /// Calculate mass matrix (linear formulation)
    /// M = integral of N^T * ρ * N dV
    virtual void calculateMassMatrix(RgMatrix& M) const override = 0;

    /// Calculate internal force vector (linear)
    /// F_int = integral of B^T * σ dV
    virtual void calculateInternalForceVector(RgVector& F) const override = 0;
};

} // namespace RgFem
