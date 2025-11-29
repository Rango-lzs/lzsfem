#pragma once

#include "RgNLSolid2dElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/**
 * @class RgNLQuad4Element
 * @brief 4-node bilinear quadrilateral element with geometric nonlinearity
 *
 * This element extends RgNLSolid2dElement with geometric nonlinearity support
 * for large deformation and large displacement analysis in 2D.
 *
 * Features:
 * - Bilinear shape functions in parametric domain r,s ∈ [-1, 1]
 * - 2×2 Gauss quadrature (4 points)
 * - Deformation gradient F computation
 * - Green-Lagrange strain formulation (material description)
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness matrix
 * - Support for hyperelastic and elastoplastic materials
 *
 * Node numbering:
 *   3 --- 2
 *   |     |
 *   0 --- 1
 */
class RgNLQuad4Element : public RgNLSolid2dElement
{
public:
    static constexpr int kNodeCount = 4;

    /// Default constructor
    RgNLQuad4Element();

    /// Parameterized constructor with node IDs
    explicit RgNLQuad4Element(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgNLQuad4Element(const RgNLQuad4Element& other);

    /// Assignment operator
    RgNLQuad4Element& operator=(const RgNLQuad4Element& other);

    /// Destructor
    virtual ~RgQuad4GeomNLElement();

    /// Return clone of this element
    virtual RgElement* clone() const override;

    /// Return the element type name
    virtual std::string typeName() const override;

    // ========== Shape Function Methods ==========
    /**
     * @brief Evaluate bilinear shape function at natural coordinates
     * @param nodeId Node index (0-3)
     * @param r Natural coordinate in x-direction [-1, 1]
     * @param s Natural coordinate in y-direction [-1, 1]
     * @param t Not used for 2D element
     * @return Shape function value
     */
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /**
     * @brief Evaluate shape function derivatives at natural coordinates
     * @param nodeId Node index (0-3)
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] dNdr ∂N/∂r
     * @param[out] dNds ∂N/∂s
     * @param[out] dNdt ∂N/∂t (always 0)
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
     * @param[out] J 3×3 Jacobian matrix
     */
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /**
     * @brief Evaluate Jacobian determinant at natural coordinates
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @return Determinant of Jacobian
     */
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /**
     * @brief Evaluate inverse of Jacobian matrix
     * @param r Natural coordinate in x-direction
     * @param s Natural coordinate in y-direction
     * @param t Not used for 2D element
     * @param[out] Jinv Inverse of Jacobian matrix
     */
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    /// Return the number of Gauss points (4 for 2×2 quadrature)
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
     * @return Area of the quadrilateral element
     */
    double getElementArea() const;

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

#endif // RGNLQUAD4ELEMENT_H
