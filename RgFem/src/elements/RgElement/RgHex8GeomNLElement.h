#ifndef RGHEX8GEOMNLELEMENT_H
#define RGHEX8GEOMNLELEMENT_H

#include "RgNLSolid3dElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/*
    RgHex8GeomNLElement
    - 8-node linear hexahedral element with geometric nonlinearity
    - Derives from RgNLSolid3dElement (nonlinear 3D solid element base)
    - Used for large deformation and large displacement analysis
    - Accounts for large deformations and rotations
    - Features:
      * Own shape function implementations (8-node linear hex)
      * 8-point Gauss quadrature (2×2×2)
      * Deformation gradient F computation
      * Green-Lagrange strain (material description)
      * Cauchy stress (spatial description)
      * Geometric (initial stress) stiffness matrix
      * Supports hyperelastic and elastoplastic materials
*/
class RgHex8GeomNLElement : public RgNLSolid3dElement
{
public:
    static constexpr int kNodeCount = 8;

    // Constructors and Destructors
    RgHex8GeomNLElement();
    explicit RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds);
    RgHex8GeomNLElement(const RgHex8GeomNLElement& other);
    RgHex8GeomNLElement& operator=(const RgHex8GeomNLElement& other);
    virtual ~RgHex8GeomNLElement();

    // Element identification
    virtual RgElement* clone() const override;
    virtual std::string typeName() const override;

    // ========== Shape Function Methods (Own Implementation) ==========
    // Shape function N_i
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;
    
    // Natural derivatives ∂N_i/∂r, ∂N_i/∂s, ∂N_i/∂t
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                   double& dNdr, double& dNds, double& dNdt) const override;

    // ========== Coordinate/Jacobian Methods (Own Implementation) ==========
    // Evaluate coordinates at natural point: x = sum(N_i * x_i)
    virtual void evaluateCoordinates(double r, double s, double t,
                                      std::array<double, 3>& coord) const override;

    // Evaluate Jacobian matrix: J[i][j] = ∂x_i/∂ξ_j
    virtual void evaluateJacobian(double r, double s, double t,
                                   std::array<std::array<double, 3>, 3>& J) const override;

    // Evaluate Jacobian determinant
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    // Evaluate inverse of Jacobian matrix
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                         std::array<std::array<double, 3>, 3>& Jinv) const override;

    // ========== Strain/Stress Calculation Methods ==========
    // Compute deformation gradient F = I + ∂u/∂X
    virtual void computeDeformationGradient(
        int gaussPointIndex,
        const std::vector<double>& displacement,
        std::array<std::array<double, 3>, 3>& F) const override;

    // Compute Green-Lagrange strain E = 0.5(C - I) where C = F^T * F
    virtual void computeGreenLagrangeStrain(
        const std::array<std::array<double, 3>, 3>& F,
        std::array<std::array<double, 3>, 3>& E) const override;

    // Compute Cauchy (true) stress from deformation and material
    virtual void computeCauchyStress(
        const std::array<std::array<double, 3>, 3>& F,
        const RgMaterial& material,
        std::array<std::array<double, 3>, 3>& sigma) const override;

    // ========== Matrix Assembly Methods ==========
    // Calculate tangent stiffness matrix (for nonlinear iteration)
    virtual void calculateTangentStiffnessMatrix(RgMatrix& Kt) const override;

    // Calculate mass matrix (nonlinear, constant in Lagrangian)
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    // Calculate internal force vector (nonlinear case)
    virtual void calculateInternalForceVector(RgVector& F) const override;

    // Calculate geometric (initial stress) stiffness
    virtual void calculateGeometricStiffnessMatrix(RgMatrix& Kg) const override;

    // ========== State Management ==========
    // Update displacement state for geometric nonlinearity tracking
    virtual void updateDisplacementState(const std::vector<double>& displacement) override;

private:
    // Hex-specific shape function helpers
    double N_linear(int nodeId, double r, double s, double t) const;
    void dN_linear(int nodeId, double r, double s, double t,
                   double& dNdr, double& dNds, double& dNdt) const;

    // Utility: matrix determinant and inverse
    static double matrixDeterminant(const std::array<std::array<double, 3>, 3>& A);
    static void matrixInverse(const std::array<std::array<double, 3>, 3>& A,
                              std::array<std::array<double, 3>, 3>& Ainv);

    // Gauss quadrature points and weights (8-point for hex)
    static constexpr int kGaussPoints = 8;
    static constexpr int kGaussPointsPerDir = 2;
    
    // Local gauss quadrature data
    static const std::array<double, 2> gaussPoints_1D;  // {-1/√3, 1/√3}
    static const std::array<double, 2> gaussWeights_1D; // {1, 1}
};

} // namespace RgFem

#endif // RGHEX8GEOMNLELEMENT_H