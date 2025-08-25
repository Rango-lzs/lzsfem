#ifndef RGNLSOLIDELEMENT_H
#define RGNLSOLIDELEMENT_H

#include "RgSolid3dElement.h"      // Base class for 3D solid elements
#include "RgFemTypedefs.h"          // Project-specific types (RgMatrix, RgVector)
#include <array>
#include <vector>
#include <memory>

namespace RgFem {

/**
 * RgNLSolidElement
 * - Base class for nonlinear 3D solid elements
 * - Handles material nonlinearity (plasticity, hyperelasticity) and/or geometric nonlinearity (large deformations)
 * - Requires material model with nonlinear stress-strain response
 */
class RgNLSolidElement : public RgSolid3dElement {
public:
    static constexpr int kMaxNodes = 8;  // Default to 8-node hex for illustration

    // Constructors
    RgNLSolidElement();
    explicit RgNLSolidElement(const std::array<int, kMaxNodes>& nodeIds);
    RgNLSolidElement(const RgNLSolidElement& other);
    virtual ~RgNLSolidElement();

    // RgElement overrides
    virtual RgElement* clone() const override = 0;  // Pure virtual for polymorphic cloning
    virtual std::string typeName() const override = 0;

    // Nonlinear mechanics interfaces
    virtual void computeStiffness(RgMatrix& Ke, const RgMaterial& material, int integrationOrder = 2) const override = 0;
    virtual void computeMass(RgMatrix& Me, double density, int integrationOrder = 2) const override;

    // Compute tangent stiffness matrix for nonlinear analysis
    virtual void computeTangentStiffness(RgMatrix& Kt, const RgMaterial& material, int integrationOrder = 2) const = 0;

    // Compute strain and stress at Gauss points (nonlinear version)
    virtual void computeStrainStressAtGauss(
        const std::array<std::array<double, 3>, kMaxNodes>& coords,
        const std::array<std::array<double, 3>, kMaxNodes>& dN_dxi,
        const std::array<double, 3 * kMaxNodes>& nodalDisp,
        const RgMaterial& material,
        std::array<double, 6>& strain,  // Strain in Voigt notation
        std::array<double, 6>& stress) const = 0;

    // Update material state (e.g., plasticity, damage)
    virtual void updateMaterialState(
        const std::array<std::array<double, 3>, kMaxNodes>& coords,
        const std::array<std::array<double, 3>, kMaxNodes>& dN_dx,
        const RgMaterial& material,
        double timeStep) const;

    // Accessors for nodal displacements
    const std::array<std::array<double, 3>, kMaxNodes>& getDisplacement() const {
        return m_currentDisplacement;
    }

    void setDisplacement(const std::array<std::array<double, 3>, kMaxNodes>& disp) {
        m_currentDisplacement = disp;
    }

protected:
    // Current nodal displacements (for state-dependent computations)
    mutable std::array<std::array<double, 3>, kMaxNodes> m_currentDisplacement;

    // Material state storage (e.g., plastic strain, damage variable)
    mutable std::vector<double> m_materialState;  // One entry per Gauss point
};

} // namespace RgFem

#endif // RGNLSOLIDELEMENT_H