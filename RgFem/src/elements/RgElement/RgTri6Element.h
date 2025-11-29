#pragma once

#include "RgLinearSolid2dElement.h"

namespace RgFem {

/**
 * @class RgTri6Element
 * @brief 6-node quadratic triangular element for 2D solid analysis
 *
 * This element uses quadratic shape functions based on barycentric (area) coordinates
 * over a triangular parametric domain. It is suitable for improved accuracy in
 * small-strain linear analysis.
 *
 * Supports:
 * - Plane stress analysis
 * - Plane strain analysis
 * - Axisymmetric analysis
 *
 * The element has 6 nodes: 3 corner nodes and 3 mid-edge nodes, each with 2 DOF (displacement in x, y).
 * Node numbering follows the standard convention:
 *        2
 *       / \
 *      5   4
 *     /     \
 *    0 - 3 - 1
 *
 * Gauss integration uses 3-point or 6-point quadrature for accurate integration of quadratic elements.
 * Default is 6-point for better accuracy.
 *
 * @author
 * @version 1.0
 */
class RgTri6Element : public RgLinearSolid2dElement
{
public:
    static constexpr int kNodeCount = 6;  ///< Number of nodes (3 corner + 3 mid-edge)

    /// Default constructor
    RgTri6Element();

    /// Parameterized constructor with node IDs
    explicit RgTri6Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgTri6Element(const RgTri6Element& other);

    /// Assignment operator
    RgTri6Element& operator=(const RgTri6Element& other);

    /// Destructor
    virtual ~RgTri6Element();

    /// Return the element type
    virtual ElementType elementType() const override;

    /// Return the element shape
    virtual ElementShape elementShape() const override;

    /// Return the element classification
    virtual ElementCategory elementClass() const override;

    /// Return the number of Gauss points (6-point for quadratic triangle)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits (Gauss points, weights, shape function derivatives)
    virtual void initTraits() override;

    /**
     * @brief Evaluate quadratic triangular shape function at natural coordinates
     * @param nodeId Node index (0-5)
     * @param r Natural coordinate (local coordinate 1)
     * @param s Natural coordinate (local coordinate 2)
     * @param t Not used for 2D element (set to 0)
     * @return Shape function value
     *
     * Corner nodes (0,1,2): Quadratic in barycentric coordinates
     * - Node 0: L1(2*L1 - 1) where L1 = 1-r-s
     * - Node 1: L2(2*L2 - 1) where L2 = r
     * - Node 2: L3(2*L3 - 1) where L3 = s
     * 
     * Mid-edge nodes (3,4,5): Linear product of barycentric coordinates
     * - Node 3: 4*L1*L2 (edge 0-1)
     * - Node 4: 4*L2*L3 (edge 1-2)
     * - Node 5: 4*L3*L1 (edge 2-0)
     */
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /**
     * @brief Evaluate derivatives of shape functions with respect to natural coordinates
     * @param nodeId Node index (0-5)
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @param[out] dNdr ∂N/∂r
     * @param[out] dNds ∂N/∂s
     * @param[out] dNdt ∂N/∂t (always 0 for 2D)
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
     * @param[out] J 3×3 Jacobian matrix (2D active submatrix is 2×2)
     */
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /**
     * @brief Evaluate Jacobian determinant at natural coordinates
     * @param r Natural coordinate
     * @param s Natural coordinate
     * @param t Not used for 2D element
     * @return Determinant of the Jacobian matrix
     */
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /**
     * @brief Compute the element area using numerical integration
     * @return Area of the triangular element
     */
    double getElementArea() const;

    /**
     * @brief Evaluate inverse of Jacobian matrix
     * @param r Natural coordinate
     * @param s Natural coordinate
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
     * @brief Helper: Get barycentric coordinate L1 = 1-r-s
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Barycentric L1
     */
    double L1(double r, double s) const { return 1.0 - r - s; }

    /**
     * @brief Helper: Get barycentric coordinate L2 = r
     * @param r Natural r-coordinate
     * @return Barycentric L2
     */
    double L2(double r) const { return r; }

    /**
     * @brief Helper: Get barycentric coordinate L3 = s
     * @param s Natural s-coordinate
     * @return Barycentric L3
     */
    double L3(double s) const { return s; }

    /**
     * @brief Helper: Get corner shape function
     * @param cornerIdx Corner node index (0-2)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Shape function value for corner node
     */
    double N_corner(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get mid-edge shape function
     * @param edgeIdx Edge index (3-5)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return Shape function value for mid-edge node
     */
    double N_midedge(int edgeIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of corner shape function with respect to r
     * @param cornerIdx Corner node index (0-2)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂r
     */
    double dN_corner_dr(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of corner shape function with respect to s
     * @param cornerIdx Corner node index (0-2)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂s
     */
    double dN_corner_ds(int cornerIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of mid-edge shape function with respect to r
     * @param edgeIdx Edge index (3-5)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂r
     */
    double dN_midedge_dr(int edgeIdx, double r, double s) const;

    /**
     * @brief Helper: Get derivative of mid-edge shape function with respect to s
     * @param edgeIdx Edge index (3-5)
     * @param r Natural r-coordinate
     * @param s Natural s-coordinate
     * @return ∂N/∂s
     */
    double dN_midedge_ds(int edgeIdx, double r, double s) const;
};

} // namespace RgFem
