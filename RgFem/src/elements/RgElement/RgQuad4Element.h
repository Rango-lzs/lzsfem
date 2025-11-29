#pragma once

#include "RgLinearSolid2dElement.h"

namespace RgFem {

/**
 * @class RgQuad4Element
 * @brief 4-node bilinear quadrilateral element for 2D solid analysis
 *
 * This element uses bilinear shape functions over the parametric domain r,s ∈ [-1, 1].
 * It is suitable for small-strain linear analysis and supports:
 * - Plane stress analysis
 * - Plane strain analysis
 * - Axisymmetric analysis
 *
 * The element has 4 corner nodes, each with 2 DOF (displacement in x, y directions).
 * Node numbering follows the standard convention:
 *   3 --- 2
 *   |     |
 *   0 --- 1
 *
 * Gauss integration uses 2×2 quadrature (4 points) for proper accuracy.
 *
 * @author
 * @version 1.0
 */
class RgQuad4Element : public RgLinearSolid2dElement
{
public:
    static constexpr int kNodeCount = 4;  ///< Number of nodes

    /// Default constructor
    RgQuad4Element();

    /// Parameterized constructor with node IDs
    explicit RgQuad4Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgQuad4Element(const RgQuad4Element& other);

    /// Assignment operator
    RgQuad4Element& operator=(const RgQuad4Element& other);

    /// Destructor
    virtual ~RgQuad4Element();

    /// Return the element type
    virtual ElementType elementType() const override;

    /// Return the element shape
    virtual ElementShape elementShape() const override;

    /// Return the element classification
    virtual ElementCategory elementClass() const override;

    /// Return the number of Gauss points (2×2 = 4 for bilinear)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits (Gauss points, weights, shape function derivatives)
    virtual void initTraits() override;

    /**
     * @brief Evaluate bilinear shape function at natural coordinates
     * @param nodeId Node index (0-3)
     * @param r Natural coordinate in x-direction [-1, 1]
     * @param s Natural coordinate in y-direction [-1, 1]
     * @param t Not used for 2D element (set to 0)
     * @return Shape function value
     */
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /**
     * @brief Evaluate derivatives of shape functions with respect to natural coordinates
     * @param nodeId Node index (0-3)
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] dNdr ∂N/∂r
     * @param[out] dNds ∂N/∂s
     * @param[out] dNdt ∂N/∂t (always 0 for 2D)
     */
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                 double& dNdr, double& dNds, double& dNdt) const override;

    /**
     * @brief Evaluate physical coordinates at natural coordinates
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] coord Physical coordinates (x, y, z)
     */
    virtual void evaluateCoordinates(double r, double s, double t,
                                    std::array<double, 3>& coord) const override;

    /**
     * @brief Evaluate Jacobian matrix at natural coordinates
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] J 3×3 Jacobian matrix (2D active submatrix is 2×2)
     */
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /**
     * @brief Evaluate Jacobian determinant at natural coordinates
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @return Determinant of the Jacobian matrix
     */
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /**
     * @brief Compute the element area
     * @return Area of the quadrilateral element
     */
    double getElementArea() const;

    /**
     * @brief Evaluate inverse of Jacobian matrix
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] Jinv Inverse of Jacobian matrix
     */
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    /**
     * @brief Calculate element stiffness matrix
     * @param[out] K Element stiffness matrix
     */
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;

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

private:
    /**
     * @brief Helper: Get bilinear shape function for node in parametric domain
     * @param nodeId Node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Shape function value
     */
    double N(int nodeId, double r, double s) const;

    /**
     * @brief Helper: Get derivative of bilinear shape function with respect to r
     * @param nodeId Node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂r
     */
    double dNdr(int nodeId, double r, double s) const;

    /**
     * @brief Helper: Get derivative of bilinear shape function with respect to s
     * @param nodeId Node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂s
     */
    double dNds(int nodeId, double r, double s) const;
};

} // namespace RgFem
