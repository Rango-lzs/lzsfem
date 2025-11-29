#pragma once

#include "RgLinearShellElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/**
 * @class RgShell4Element
 * @brief 4-node bilinear shell element (quadrilateral)
 *
 * This element combines in-plane bilinear behavior with plate bending (Kirchhoff shell theory).
 *
 * Features:
 * - Bilinear shape functions in-plane (r,s ∈ [-1, 1])
 * - Cubic shape functions for transverse shear (Mindlin plate theory optional)
 * - 6 DOF per node (3 translations + 3 rotations)
 * - 2×2 Gauss quadrature in-plane
 * - Multiple through-thickness integration points
 * - Membrane + bending stiffness
 * - Shear stiffness
 * - Support for laminated composite shells
 *
 * Node numbering:
 *   3 --- 2
 *   |     |
 *   0 --- 1
 *
 * Local coordinate system: x-y plane for in-plane, z for thickness
 */
class RgShell4Element : public RgLinearShellElement
{
public:
    static constexpr int kNodeCount = 4;
    static constexpr int kDofsPerNode = 6;  // 3 translations + 3 rotations

    /// Default constructor
    RgShell4Element();

    /// Parameterized constructor with node IDs
    explicit RgShell4Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgShell4Element(const RgShell4Element& other);

    /// Assignment operator
    RgShell4Element& operator=(const RgShell4Element& other);

    /// Destructor
    virtual ~RgShell4Element();

    /// Return clone of this element
    virtual RgElement* clone() const override;

    /// Return the element type name
    virtual std::string typeName() const override;

    // ========== Shape Function Methods ==========
    /// Evaluate bilinear shape function at natural coordinates
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

    /// Return the number of in-plane Gauss points (4 for 2×2)
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

    /// Get the element area (projected on xy-plane)
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

    /// Helper: Get bilinear shape function
    double N(int nodeId, double r, double s) const;

    /// Helper: Get derivative of bilinear shape function with respect to r
    double dNdr(int nodeId, double r, double s) const;

    /// Helper: Get derivative of bilinear shape function with respect to s
    double dNds(int nodeId, double r, double s) const;
};

} // namespace RgFem

#endif // RGSHELL4ELEMENT_H
