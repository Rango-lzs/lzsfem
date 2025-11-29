#pragma once

#include "RgLinearShellElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/**
 * @class RgShell3Element
 * @brief 3-node triangular shell element
 *
 * This element combines in-plane linear behavior with plate bending (Kirchhoff shell theory).
 *
 * Features:
 * - Linear shape functions in-plane (barycentric coordinates)
 * - Cubic shape functions for transverse shear (Mindlin plate theory optional)
 * - 6 DOF per node (3 translations + 3 rotations)
 * - 1-point Gauss quadrature in-plane (constant strain triangle)
 * - Multiple through-thickness integration points
 * - Membrane + bending stiffness
 * - Shear stiffness
 * - Support for laminated composite shells
 *
 * Node numbering:
 *      2
 *     / \
 *    0 - 1
 *
 * Local coordinate system: x-y plane for in-plane, z for thickness
 */
class RgShell3Element : public RgLinearShellElement
{
public:
    static constexpr int kNodeCount = 3;
    static constexpr int kDofsPerNode = 6;  // 3 translations + 3 rotations

    /// Default constructor
    RgShell3Element();

    /// Parameterized constructor with node IDs
    explicit RgShell3Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgShell3Element(const RgShell3Element& other);

    /// Assignment operator
    RgShell3Element& operator=(const RgShell3Element& other);

    /// Destructor
    virtual ~RgShell3Element();

    /// Return clone of this element
    virtual RgElement* clone() const override;

    /// Return the element type name
    virtual std::string typeName() const override;

    // ========== Shape Function Methods ==========
    /// Evaluate linear triangular shape function at natural coordinates
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /// Evaluate shape function derivatives at natural coordinates
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                 double& dNdr, double& dNds, double& dNdt) const override;

    /// Evaluate physical coordinates at natural coordinates
    virtual void evaluateCoordinates(double r, double s, double t,
                                    std::array<double, 3>& coord) const override;

    /// Evaluate Jacobian matrix at natural coordinates
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /// Evaluate Jacobian determinant at natural coordinates
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /// Evaluate inverse of Jacobian matrix
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    /// Return the number of Gauss points (1 for constant strain triangle)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits
    virtual void initTraits() override;

    // ========== Shell-Specific Methods ==========
    /// Get the shell normal vector at a point
    virtual void getShellNormal(double r, double s, std::array<double, 3>& normal) const;

    /// Get the shell thickness
    virtual double getShellThickness() const;

    /// Set the shell thickness
    virtual void setShellThickness(double thickness);

    /// Get the element area (projected on surface)
    virtual double getElementArea() const;

    // ========== Matrix Assembly Methods ==========
    /// Calculate element stiffness matrix
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;

    /// Calculate element mass matrix
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    /// Calculate internal force vector
    virtual void calculateInternalForceVector(RgVector& F) const override;

private:
    double m_thickness = 1.0;  ///< Shell thickness

    /// Helper: Get barycentric coordinate L1 = 1-r-s
    double L1(double r, double s) const { return 1.0 - r - s; }

    /// Helper: Get barycentric coordinate L2 = r
    double L2(double r) const { return r; }

    /// Helper: Get barycentric coordinate L3 = s
    double L3(double s) const { return s; }
};

} // namespace RgFem

#endif // RGSHELL3ELEMENT_H
