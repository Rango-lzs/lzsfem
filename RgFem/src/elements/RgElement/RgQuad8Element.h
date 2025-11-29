#pragma once

#include "RgLinearSolid2dElement.h"

namespace RgFem {

/**
 * @class RgQuad8Element
 * @brief 8-node serendipity quadrilateral element for 2D solid analysis
 *
 * This element uses quadratic serendipity shape functions (mid-side nodes but no center node)
 * over the parametric domain r,s ∈ [-1, 1]. It is suitable for improved accuracy in
 * small-strain linear analysis while maintaining computational efficiency.
 *
 * Supports:
 * - Plane stress analysis
 * - Plane strain analysis
 * - Axisymmetric analysis
 *
 * The element has 8 nodes: 4 corner nodes and 4 mid-edge nodes, each with 2 DOF (displacement in x, y).
 * Node numbering follows the standard convention:
 *   3 --- 6 --- 2
 *   |           |
 *   7     *     5
 *   |           |
 *   0 --- 4 --- 1
 *
 * Gauss integration uses 3×3 quadrature (9 points) for accurate integration.
 *
 * @author
 * @version 1.0
 */
class RgQuad8Element : public RgLinearSolid2dElement
{
public:
    static constexpr int kNodeCount = 8;  ///< Number of nodes (4 corner + 4 mid-edge)

    /// Default constructor
    RgQuad8Element();

    /// Parameterized constructor with node IDs
    explicit RgQuad8Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgQuad8Element(const RgQuad8Element& other);

    /// Assignment operator
    RgQuad8Element& operator=(const RgQuad8Element& other);

    /// Destructor
    virtual ~RgQuad8Element();

    /// Return the element type
    virtual ElementType elementType() const override;

    /// Return the element shape
    virtual ElementShape elementShape() const override;

    /// Return the element classification
    virtual ElementCategory elementClass() const override;

    /// Return the number of Gauss points (3×3 = 9 for quadratic)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits (Gauss points, weights, shape function derivatives)
    virtual void initTraits() override;

    /**
     * @brief Evaluate serendipity quadratic shape function at natural coordinates
     * @param nodeId Node index (0-7)
     * @param r Natural coordinate in x-direction [-1, 1]
     * @param s Natural coordinate in y-direction [-1, 1]
     * @param t Not used for 2D element (set to 0)
     * @return Shape function value
     *
     * Corner nodes (0,1,2,3): quadratic in r and s, proportional to (r±1)(s±1)(r/s-1/2)
     * Mid-edge nodes (4,5,6,7): quadratic in one direction, linear in other
     */
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /**
     * @brief Evaluate derivatives of shape functions with respect to natural coordinates
     * @param nodeId Node index (0-7)
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
     * @brief Compute the element area using numerical integration
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
     * @brief Helper: Get serendipity quadratic shape function for corner node
     * @param cornerIdx Corner node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Shape function value
     */
    double N_corner(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get serendipity quadratic shape function for mid-edge node
     * @param edgeIdx Edge index (4-7 for mid-edge nodes)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Shape function value
     */
    double N_midedge(int edgeIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of corner shape function with respect to r
     * @param cornerIdx Corner node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂r
     */
    double dN_corner_dr(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of corner shape function with respect to s
     * @param cornerIdx Corner node index (0-3)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂s
     */
    double dN_corner_ds(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of mid-edge shape function with respect to r
     * @param edgeIdx Edge index (4-7)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂r
     */
    double dN_midedge_dr(int edgeIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of mid-edge shape function with respect to s
     * @param edgeIdx Edge index (4-7)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂s
     */
    double dN_midedge_ds(int edgeIdx, double r, double s) const;
};

} // namespace RgFem
