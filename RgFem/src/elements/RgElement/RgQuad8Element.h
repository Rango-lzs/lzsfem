#pragma once

#include "RgLinearSolid2dElement.h"



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
    virtual ElementCategory elementCategory() const override;

    /**
     * @brief Calculate element stiffness matrix
     * @param[out] K Element stiffness matrix
     */
    virtual void calculateStiffnessMatrix(Matrix& K) override;

    /**
     * @brief Calculate element mass matrix
     * @param[out] M Element mass matrix
     */
    virtual void calculateMassMatrix(Matrix& M) override;

    /**
     * @brief Calculate internal force vector
     * @param[out] F Element internal force vector
     */
    virtual void calculateInternalForceVector(RgVector& F) override;
};