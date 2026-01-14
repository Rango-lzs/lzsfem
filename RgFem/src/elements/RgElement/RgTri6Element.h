#pragma once

#include "RgLinearSolid2dElement.h"



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
    virtual ElementCategory elementCategory() const override;


    /**
     * @brief Calculate element stiffness matrix
     * @param[out] K Element stiffness matrix
     */
    virtual void calculateStiffnessMatrix(RgMatrix& K) override;

    /**
     * @brief Calculate element mass matrix
     * @param[out] M Element mass matrix
     */
    virtual void calculateMassMatrix(RgMatrix& M) override;

    /**
     * @brief Calculate internal force vector
     * @param[out] F Element internal force vector
     */
    virtual void calculateInternalForceVector(RgVector& F)  override;

};