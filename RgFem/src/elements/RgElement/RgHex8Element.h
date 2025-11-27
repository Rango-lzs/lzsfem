#ifndef RGHEX8ELEMENT_H
#define RGHEX8ELEMENT_H

#include "RgSolid3dElement.h"
#include <array>
#include <vector>
#include <string>

namespace RgFem {

// forward declarations for project types
class RgMatrix;
class RgVector;
class RgMaterial;
class FEMaterialPoint;
class DumpStream;
class Matrix3ds;

/*
    RgHex8Element
    - 8-node linear hexahedral solid element (Q1 element)
    - Isoparametric element with trilinear shape functions
    - Derivatives from RgSolid3dElement
    - Features:
      * Reduced integration (1 Gauss point) or full integration (8 Gauss points)
      * Supports small strain and geometrically nonlinear analysis
      * Suitable for 3D continuum mechanics problems
*/
class RgHex8Element : public RgSolid3dElement {
public:
    static constexpr int kNodeCount = 8;
    static constexpr int kNumFaces = 6;
    static constexpr int kNumEdges = 12;

    // Constructors and Destructors
    RgHex8Element();
    explicit RgHex8Element(const std::array<int, kNodeCount>& nodeIds);
    RgHex8Element(const RgHex8Element& other);
    RgHex8Element& operator=(const RgHex8Element& other);
    virtual ~RgHex8Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementClass() const;
    
    // Element Properties
    virtual int getNumberOfNodes() const override { return kNodeCount; }
    virtual int getNumberOfGaussPoints() const;
    virtual int getNumberOfFaces() const { return kNumFaces; }
    virtual int getNumberOfEdges() const { return kNumEdges; }

    // Initialize element traits (shape functions, integration rules)
    virtual void initTraits() override;

    // Geometric operations
    virtual Vector3d evaluateCoordinates(const Vector3d& naturalCoord) const;
    virtual Matrix3d evaluateJacobian(const Vector3d& naturalCoord) const;
    double evaluateJacobianDeterminant(const Vector3d& naturalCoord) const;
    virtual Matrix3d evaluateJacobianInverse(const Vector3d& naturalCoord) const;

    // Shape function evaluations
    // Evaluate shape functions at natural coordinates (r, s, t) in [-1, 1]^3
    void evaluateShapeFunctions(double r, double s, double t, std::vector<double>& N) const;
    
    // Evaluate shape function derivatives with respect to natural coordinates
    void evaluateShapeDerivatives(double r, double s, double t,
                                   std::vector<double>& dN_dr,
                                   std::vector<double>& dN_ds,
                                   std::vector<double>& dN_dt) const;
    
    // Evaluate shape function second derivatives
    void evaluateShapeDerivatives2(double r, double s, double t,
                                    std::vector<double>& d2N_drr,
                                    std::vector<double>& d2N_dss,
                                    std::vector<double>& d2N_dtt,
                                    std::vector<double>& d2N_drs,
                                    std::vector<double>& d2N_dst,
                                    std::vector<double>& d2N_drt) const;

    // Physical field evaluations at arbitrary point
    Vector3d evaluateField(const Vector3d* nodeValues, const Vector3d& naturalCoord) const;
    double evaluateScalarField(const double* nodeValues, const Vector3d& naturalCoord) const;

    // Stiffness and Mass Matrices
    virtual void calculateStiffnessMatrix(Matrix& K) const override;
    virtual void calculateMassMatrix(Matrix& M) const override;
    virtual void calculateDampingMatrix(Matrix& C) const override;
    
    // Tangent stiffness for nonlinear analysis
    virtual void calculateTangentStiffnessMatrix(Matrix& Kt) const;

    // Internal force vector (residual)
    virtual void calculateInternalForceVector(Vector& F) const;

    // Strain and Stress Calculations
    virtual void calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain) override;
    
    // Strain energy calculation
    virtual double calculateStrainEnergy() const;
    
    // Kinetic energy calculation
    virtual double calculateKineticEnergy() const;

    // Face and Edge Operations
    virtual void getFaceNodeIds(int faceId, std::array<int, 4>& faceNodes) const;
    virtual void getEdgeNodeIds(int edgeId, std::array<int, 2>& edgeNodes) const;

    // Loading Operations
    virtual void applyBodyForce(const Vector3d& force, Vector& F) const;
    virtual void applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const;
    virtual void applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const;

    // Material point access
    virtual FEMaterialPoint* getMaterialPoint(int gaussPtId);
    virtual const FEMaterialPoint* getMaterialPoint(int gaussPtId) const;

    // State management
    virtual void updateState(double timeStep);
    virtual void commitState();
    virtual void resetState();

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

    // Utility functions
    bool isValidNaturalCoordinate(const Vector3d& naturalCoord) const;
    double getCharacteristicLength() const;

protected:
    // Initialize Gauss quadrature points and weights
    void initializeGaussPoints();
    
    // Compute B matrix (strain-displacement matrix) at natural coordinates
    void computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const;
    
    // Compute geometric stiffness matrix for large deformations
    void computeGeometricStiffness(Matrix& Kg) const;

private:
    // Gauss point data
    std::vector<double> m_gaussR, m_gaussS, m_gaussT, m_gaussW;
    
    // Cached Jacobian data
    mutable std::vector<Matrix3d> m_jacobianInverse;
    mutable std::vector<double> m_jacobianDet;
};

} // namespace RgFem

#endif // RGHEX8ELEMENT_H