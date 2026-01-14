#pragma once

#include "RgNLBeamElement.h"
#include "RgFemTypedefs.h"
#include <array>
#include <vector>



// Forward declarations
class RgMaterial;

/**
 * @class RgBeam3dGeomNLElement
 * @brief 2-node 3D Timoshenko beam element with geometric nonlinearity
 *
 * This element extends RgNLBeamElement with geometric nonlinearity support
 * for large deformation and large displacement analysis in 3D.
 *
 * Features:
 * - Cubic shape functions for displacement
 * - Quadratic shape functions for rotation (through derivatives)
 * - 2-point Gauss quadrature along the beam axis
 * - Timoshenko beam theory (includes shear deformation)
 * - 6 DOF per node (3 translations + 3 rotations)
 * - Local-to-global transformation matrices
 * - Deformation gradient F computation
 * - Green-Lagrange strain formulation
 * - Cauchy stress (spatial description)
 * - Geometric (initial stress) stiffness matrix
 * - Support for hyperelastic and elastoplastic materials
 *
 * Node numbering:
 *   0 -------- 1
 *
 * Local coordinate system: x-axis along beam, y and z perpendicular
 * 
 * DOF ordering per node: [ux, uy, uz, rx, ry, rz]
 * Total DOF per element: 12 (6 per node × 2 nodes)
 */
class RgBeam3dGeomNLElement : public RgNLBeamElement
{
public:
    static constexpr int kNodeCount = 2;
    static constexpr int kDofsPerNode = 6;  // 3 translations + 3 rotations

    /// Default constructor
    RgBeam3dGeomNLElement();

    /// Parameterized constructor with node IDs
    explicit RgBeam3dGeomNLElement(const std::array<int, kNodeCount>& nodeIds);

    /// Copy constructor
    RgBeam3dGeomNLElement(const RgBeam3dGeomNLElement& other);

    /// Assignment operator
    RgBeam3dGeomNLElement& operator=(const RgBeam3dGeomNLElement& other);

    /// Destructor
    virtual ~RgBeam3dGeomNLElement();

    /// Return clone of this element
    virtual RgElement* clone() const override;

    /// Return the element type name
    virtual std::string typeName() const override;

    // ========== Shape Function Methods ==========
    /// Evaluate shape function at natural coordinate (along beam axis)
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;

    /// Evaluate shape function derivative with respect to natural coordinate
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                 double& dNdr, double& dNds, double& dNdt) const override;

    /// Evaluate physical coordinates at natural coordinate
    virtual void evaluateCoordinates(double r, double s, double t,
                                    std::array<double, 3>& coord) const override;

    /// Evaluate Jacobian matrix at natural coordinate
    virtual void evaluateJacobian(double r, double s, double t,
                                 std::array<std::array<double, 3>, 3>& J) const override;

    /// Evaluate Jacobian determinant at natural coordinate
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;

    /// Evaluate inverse of Jacobian matrix
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    /// Return the number of Gauss points (2-point for 3D beam)
    virtual int getNumberOfGaussPoints() const override;

    /// Initialize element traits
    virtual void initTraits() override;

    // ========== Beam-Specific Methods ==========
    /// Get the length of the beam
    virtual double getBeamLength() const override;

    /// Get the local-to-global transformation matrix
    virtual void getLocalToGlobalMatrix(std::array<std::array<double, 3>, 3>& T) const override;

    /// Get the local x-axis direction (along beam)
    virtual void getLocalXAxis(std::array<double, 3>& xAxis) const override;

    /// Get the local y-axis direction (perpendicular to beam)
    virtual void getLocalYAxis(std::array<double, 3>& yAxis) const override;

    /// Get the local z-axis direction (perpendicular to beam)
    virtual void getLocalZAxis(std::array<double, 3>& zAxis) const override;

    // ========== Nonlinear Analysis Methods ==========
    /// Compute deformation gradient F at a Gauss point
    virtual void computeDeformationGradient(int gaussPointIndex,
                                           const RgVector& displacement,
                                           std::array<std::array<double, 3>, 3>& F) const override;

    /// Compute displacement gradient ∇u at a Gauss point
    virtual void computeDisplacementGradient(int gaussPointIndex,
                                            const RgVector& displacement,
                                            std::array<std::array<double, 3>, 3>& dispGrad) const override;

    // ========== Matrix Assembly Methods ==========
    /// Calculate element tangent stiffness matrix
    virtual void calculateTangentStiffnessMatrix(RgMatrix& K) const override;

    /// Calculate element geometric stiffness matrix
    virtual void calculateGeometricStiffnessMatrix(RgMatrix& Kg) const override;

    /// Calculate element mass matrix
    virtual void calculateMassMatrix(RgMatrix& M) const override;

    /// Calculate internal force vector
    virtual void calculateInternalForceVector(RgVector& F) const override;

private:
    /// Helper: Cubic Hermite shape function for displacement
    double N_cubic(int nodeId, double r) const;

    /// Helper: Derivative of cubic Hermite shape function
    double dN_cubic_dr(int nodeId, double r) const;

    /// Compute local coordinate system axes
    void computeLocalAxes(std::array<double, 3>& xAxis,
                         std::array<double, 3>& yAxis,
                         std::array<double, 3>& zAxis) const;
};



#endif // RGBEAM3DGEOMNLELEMENT_H
