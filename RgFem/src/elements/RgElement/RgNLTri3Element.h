#pragma once

#include "RgNLSolid2dElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/**
 * @class RgNLTri3Element
 * @brief 3-node linear triangular element with geometric nonlinearity
 *
 * This element extends RgNLSolid2dElement with geometric nonlinearity support
 * for large deformation and large displacement analysis in 2D.
 *
 * Features:
 * - Linear triangular shape functions in barycentric coordinates
 * - 1-point Gauss quadrature (constant strain triangle)
 * - Deformation gradient F computation
 * - Green-Lagrange strain formulation (material description)
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness matrix
 * - Support for hyperelastic and elastoplastic materials
 *
 * Node numbering:
 *   2
 *  / \
 * 0 - 1
 */
class RgNLTri3Element : public RgNLSolid2dElement
{
public:
    static constexpr int kNodeCount = 3;

    /// Default constructor
    RgNLTri3Element();

    /// Parameterized constructor with node IDs
    explicit RgNLTri3Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgNLTri3Element(const RgNLTri3Element& other);

    /// Assignment operator
    RgNLTri3Element& operator=(const RgNLTri3Element& other);

    /// Destructor
    virtual ~RgTri3GeomNLElement();

    /// Return clone of this element
    virtual RgElement* clone() const override;

    /// Return the element type name
    virtual std::string typeName() const override;

    // ========== Shape Function Methods ==========
    /**
     * @brief Evaluate shape function at natural coordinates
     * @param nodeId Node index (0-2)
     * @param r Natural coordinate (local coord 1)
     * @param s Natural coordinate (local coord 2)
     * @param t Not used for 2D element
     * @return Shape function value
     */
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /**
     * @brief Evaluate shape function derivatives at natural coordinates
     * @param nodeId Node index (0-2)
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @param[out] dNdr ∂N/∂r
     * @param[out] dNds ∂N/∂s
     * @param[out] dNdt ∂N/∂t (always 0)
     */
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                 double& dNdr, double& dNds, double& dNdt) const override;

    /**
     * @brief Evaluate physical coordinates at natural coordinates
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @param[out] coord Physical coordinates (x, y, z)
     */
    virtual void evaluateCoordinates(double r, double s, double t,
                                    std::array<double, 3>& coord) const override;

    /**
     * @brief Evaluate Jacobian matrix at natural coordinates
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @param[out] J 3×3 Jacobian matrix
     */
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /**
     * @brief Evaluate Jacobian determinant at natural coordinates
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @return Determinant of Jacobian
     */
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /**
     * @brief Evaluate inverse of Jacobian matrix
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @param[out] Jinv Inverse of Jacobian matrix
     */
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    /// Return the number of Gauss points (1 for constant strain)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits
    virtual void initTraits() override;

    // ========== Nonlinear Analysis Methods ==========
    /**
     * @brief Compute deformation gradient F at a Gauss point
     * @param gaussPointIndex Index of the Gauss point
     * @param displacement Node displacement increments
     * @param[out] F 3×3 deformation gradient matrix
     */
    virtual void computeDeformationGradient(int gaussPointIndex,
                                           const RgVector& displacement,
                                           std::array<std::array<double, 3>, 3>& F) const override;

    /**
     * @brief Compute displacement gradient ∇u at a Gauss point
     * @param gaussPointIndex Index of the Gauss point
     * @param displacement Node displacement increments
     * @param[out] dispGrad 3×3 displacement gradient matrix
     */
    virtual void computeDisplacementGradient(int gaussPointIndex,
                                            const RgVector& displacement,
                                            std::array<std::array<double, 3>, 3>& dispGrad) const override;

    // ========== Matrix Assembly Methods ==========
    /**
     * @brief Calculate element tangent stiffness matrix (linear + geometric)
     * @param[out] K Element tangent stiffness matrix
     */
    virtual void calculateTangentStiffnessMatrix(RgMatrix& K) const override;

    /**
     * @brief Calculate element geometric (initial stress) stiffness matrix
     * @param[out] Kg Element geometric stiffness matrix
     */
    virtual void calculateGeometricStiffnessMatrix(RgMatrix& Kg) const override;

    /**
     * @brief Calculate element mass matrix
     * @param[out] M Element mass matrix
     */
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    /**
     * @brief Calculate internal force vector
     * @param[out] F Element internal force vector
     */
    virtual void calculateInternalForceVector(RgVector& F) const override;

    /**
     * @brief Get element area
     * @return Area of the triangular element
     */
    double getElementArea() const;

private:
    /**
     * @brief Helper: Get barycentric coordinate L1 = 1-r-s
     */
    double L1(double r, double s) const { return 1.0 - r - s; }

    /**
     * @brief Helper: Get barycentric coordinate L2 = r
     */
    double L2(double r) const { return r; }

    /**
     * @brief Helper: Get barycentric coordinate L3 = s
     */
    double L3(double s) const { return s; }
};

} // namespace RgFem

#endif // RGNLTRI3ELEMENT_H
