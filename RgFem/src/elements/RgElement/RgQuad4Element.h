#pragma once

#include "RgLinearSolid2dElement.h"



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
    virtual void calculateInternalForceVector(std::vector<double>& F) override;

    // Additional methods required by the updated base class
    virtual void calculateDampingMatrix(Matrix& C) override;
    virtual void getStress(RgMaterialPoint& matPt, StressTensor& stress) override;
    virtual void getStrain(RgMaterialPoint& matPt, StrainTensor& strain) override;
    //virtual void updateStress(RgMaterialPoint& matPt) override;
    //virtual void updateStrain(RgMaterialPoint& matPt) override;
    virtual double calculateStrainEnergy() override;
    virtual double calculateKineticEnergy() override;
};