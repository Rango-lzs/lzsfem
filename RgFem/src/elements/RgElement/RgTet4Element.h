#ifndef RGTET4ELEMENT_H
#define RGTET4ELEMENT_H

#include "RgSolid3dElement.h"
#include <array>
#include <vector>

namespace RgFem {

// forward declarations
class Matrix3ds;

/*
    RgTet4Element
    - 4-node linear tetrahedral solid element
    - Isoparametric element with linear shape functions
    - Derives from RgSolid3dElement
    - Features:
      * Constant strain element (C0 element)
      * Full integration (1 Gauss point)
      * Supports small strain analysis
      * Suitable for mesh generation and computational efficiency
*/
class RgTet4Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 4;
    static constexpr int kNumFaces = 4;
    static constexpr int kNumEdges = 6;

    // Constructors and Destructors
    RgTet4Element();
    explicit RgTet4Element(const std::array<int, kNodeCount>& nodeIds);
    RgTet4Element(const RgTet4Element& other);
    RgTet4Element& operator=(const RgTet4Element& other);
    virtual ~RgTet4Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementClass() const;

    // Element Properties
    virtual int getNumberOfNodes() const override { return kNodeCount; }
    virtual int getNumberOfGaussPoints() const;
    virtual int getNumberOfFaces() const { return kNumFaces; }
    virtual int getNumberOfEdges() const { return kNumEdges; }

    // Initialize element traits
    virtual void initTraits() override;

    // Geometric operations
    virtual Vector3d evaluateCoordinates(const Vector3d& naturalCoord) const;
    virtual Matrix3d evaluateJacobian(const Vector3d& naturalCoord) const;
    double evaluateJacobianDeterminant(const Vector3d& naturalCoord) const;
    virtual Matrix3d evaluateJacobianInverse(const Vector3d& naturalCoord) const;

    // Shape function evaluations
    void evaluateShapeFunctions(double r, double s, double t, std::vector<double>& N) const;
    void evaluateShapeDerivatives(double r, double s, double t,
                                   std::vector<double>& dN_dr,
                                   std::vector<double>& dN_ds,
                                   std::vector<double>& dN_dt) const;

    // Physical field evaluations
    Vector3d evaluateField(const Vector3d* nodeValues, const Vector3d& naturalCoord) const;
    double evaluateScalarField(const double* nodeValues, const Vector3d& naturalCoord) const;

    // FEM Matrix Calculations
    virtual void calculateStiffnessMatrix(Matrix& K) const override;
    virtual void calculateMassMatrix(Matrix& M) const override;
    virtual void calculateDampingMatrix(Matrix& C) const override;

    // Strain and Stress Calculations
    virtual void calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain) override;

    // Face and Edge Operations
    virtual void getFaceNodeIds(int faceId, std::array<int, 3>& faceNodes) const;
    virtual void getEdgeNodeIds(int edgeId, std::array<int, 2>& edgeNodes) const;

    // Loading Operations
    virtual void applyBodyForce(const Vector3d& force, Vector& F) const;
    virtual void applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const;
    virtual void applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const;

    // Material point access
    virtual FEMaterialPoint* getMaterialPoint(int gaussPtId);
    virtual const FEMaterialPoint* getMaterialPoint(int gaussPtId) const;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

    // Utility functions
    bool isValidNaturalCoordinate(const Vector3d& naturalCoord) const;
    double getCharacteristicLength() const;
    double getVolume() const;

protected:
    // Compute B matrix (strain-displacement matrix)
    void computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const;

private:
    // Gauss point data (1 point at centroid for linear tet)
    std::vector<double> m_gaussR, m_gaussS, m_gaussT, m_gaussW;
    
    // Cached Jacobian data
    mutable std::vector<Matrix3d> m_jacobianInverse;
    mutable std::vector<double> m_jacobianDet;
};

} // namespace RgFem

#endif // RGTET4ELEMENT_H
