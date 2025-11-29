#pragma once

#include "RgNLBeamElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>

namespace RgFem {

// Forward declarations
class RgMaterial;

/**
 * @class RgBeam2dGeomNLElement
 * @brief 2-node 2D Timoshenko beam element with geometric nonlinearity
 *
 * This element extends RgNLBeamElement with geometric nonlinearity support
 * for large deformation and large displacement analysis in 2D (planar).
 *
 * Features:
 * - Linear shape functions for displacement (along beam axis)
 * - Quadratic shape functions for rotation (through derivatives)
 * - 2-point Gauss quadrature along the beam axis
 * - Timoshenko beam theory (includes shear deformation)
 * - 3 DOF per node (2 translations + 1 rotation in xy-plane)
 * - Local-to-global transformation matrices
 * - Deformation gradient F computation (2D)
 * - Green-Lagrange strain formulation
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness matrix
 * - Support for hyperelastic and elastoplastic materials
 *
 * Node numbering:
 *   0 -------- 1
 *
 * Local coordinate system: x-axis along beam (in xy-plane), y perpendicular
 * 
 * DOF ordering per node: [ux, uy, rz]
 * Total DOF per element: 6 (3 per node × 2 nodes)
 *
 * Theory:
 * - Planar Timoshenko beam with shear deformation
 * - Updated Lagrangian description with geometric nonlinearity
 * - Deformation gradient: F = I + ∇u (2D: 2×2 matrix)
 * - Green-Lagrange strain: E = 0.5(F^T*F - I)
 * - Material hyperelasticity with isotropic response
 */
class RgBeam2dGeomNLElement : public RgNLBeamElement
{
public:
    static constexpr int kNodeCount = 2;
    static constexpr int kDofsPerNode = 3;  // ux, uy, rz
    static constexpr int kTotalDofs = kNodeCount * kDofsPerNode;

    /// Default constructor
    RgBeam2dGeomNLElement();

    /// Parameterized constructor with node IDs
    explicit RgBeam2dGeomNLElement(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgBeam2dGeomNLElement(const RgBeam2dGeomNLElement& other);

    /// Assignment operator
    RgBeam2dGeomNLElement& operator=(const RgBeam2dGeomNLElement& other);

    /// Destructor
    virtual ~RgBeam2dGeomNLElement();

    /// Clone this element
    virtual RgElement* clone() const override;

    /// Return element type name
    virtual std::string typeName() const override;

    /// Get number of nodes
    virtual int getNumberOfNodes() const override { return kNodeCount; }

    /// Get number of Gauss points (2-point along beam)
    virtual int getNumberOfGaussPoints() const override;

    /// Get number of DOFs per node
    virtual int getDofsPerNode() const override { return kDofsPerNode; }

    /// Get total DOFs in element
    virtual int getTotalDofs() const override { return kTotalDofs; }

    /// Evaluate beam length
    virtual double evaluateLength() const override;

    /// Get beam axis vector
    virtual Vector3d getAxis() const override;

    // Shape function evaluations at parametric coordinate r ∈ [-1, 1]
    /// Evaluate shape functions N(r) = [N0(r), N1(r)]
    virtual void evaluateShapeFunctions(double r, std::vector<double>& N) const override;

    /// Evaluate shape function derivatives dN/dr
    virtual void evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const override;

    // Jacobian operations
    /// Evaluate Jacobian matrix at parameter r
    void evaluateJacobian(double r, std::array<std::array<double, 2>, 2>& J) const;

    /// Evaluate Jacobian determinant
    double evaluateJacobianDeterminant(double r) const;

    /// Evaluate Jacobian inverse
    void evaluateJacobianInverse(double r, std::array<std::array<double, 2>, 2>& Jinv) const;

    // Coordinate evaluation
    /// Evaluate physical coordinates at parameter r
    void evaluateCoordinates(double r, std::array<double, 3>& coord) const;

    // Local coordinate system
    /// Compute local beam axes (x along beam, y perpendicular in xy-plane)
    void computeLocalAxes(std::array<double, 2>& localX, std::array<double, 2>& localY) const;

    /// Get local-to-global rotation matrix (2×2)
    void getLocalToGlobalMatrix(std::array<std::array<double, 2>, 2>& R) const;

    // Nonlinear analysis methods
    /// Compute deformation gradient F at a Gauss point (2D: 2×2 matrix)
    /// @param gaussPointIndex Index of the Gauss point
    /// @param displacement Nodal displacement vector [ux0, uy0, rz0, ux1, uy1, rz1]
    /// @param F Output 2×2 deformation gradient matrix
    virtual void computeDeformationGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                            std::array<std::array<double, 3>, 3>& F) const override;

    /// Compute displacement gradient ∇u at a Gauss point (2D: 2×2 matrix)
    /// @param gaussPointIndex Index of the Gauss point
    /// @param displacement Nodal displacement vector [ux0, uy0, rz0, ux1, uy1, rz1]
    /// @param dispGrad Output 2×2 displacement gradient matrix
    virtual void computeDisplacementGradient(int gaussPointIndex, const std::vector<double>& displacement,
                                             std::array<std::array<double, 3>, 3>& dispGrad) const override;

    // FEM matrix calculations
    /// Calculate stiffness matrix K (nonlinear formulation with geometric stiffness)
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;

    /// Calculate mass matrix M
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    /// Calculate internal force vector F_int
    virtual void calculateInternalForceVector(RgVector& F) const override;

private:
    /// Helper function: Linear shape function N0(r) = (1-r)/2
    double N_linear(int nodeId, double r) const;

    /// Helper function: Derivative of linear shape function dN/dr
    double dN_linear_dr(int nodeId) const;

    /// Helper function: Get Gauss point coordinates and weights
    void getGaussPointData(int gaussPointIndex, double& r, double& weight) const;
};

} // namespace RgFem
