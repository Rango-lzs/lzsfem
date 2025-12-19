#ifndef RGHEX20ELEMENT_H
#define RGHEX20ELEMENT_H

#include "RgSolid3dElement.h"
#include <array>
#include <vector>

namespace RgFem {

// forward declarations
class Matrix3ds;

/*
    RgHex20Element
    - 20-node serendipity hexahedral solid element (Q2 element)
    - Isoparametric element with quadratic shape functions
    - Derives from RgSolid3dElement
    - Features:
      * Biquadratic-bilinear interpolation
      * Full integration (8 Gauss points)
      * Supports small strain analysis
      * Better accuracy for curved boundaries
      * Node distribution: 8 corner + 12 mid-edge nodes
*/
class RgHex20Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 20;
    static constexpr int kNumFaces = 6;
    static constexpr int kNumEdges = 12;

    // Constructors and Destructors
    RgHex20Element();
    explicit RgHex20Element(const std::array<int, kNodeCount>& nodeIds);
    RgHex20Element(const RgHex20Element& other);
    RgHex20Element& operator=(const RgHex20Element& other);
    virtual ~RgHex20Element();

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
    void evaluateShapeDerivatives2(double r, double s, double t,
                                    std::vector<double>& d2N_drr,
                                    std::vector<double>& d2N_dss,
                                    std::vector<double>& d2N_dtt,
                                    std::vector<double>& d2N_drs,
                                    std::vector<double>& d2N_dst,
                                    std::vector<double>& d2N_drt) const;

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

    // Face and Edge Operations (8-node quad faces, 3-node edges)
    virtual void getFaceNodeIds(int faceId, std::array<int, 8>& faceNodes) const;
    virtual void getEdgeNodeIds(int edgeId, std::array<int, 3>& edgeNodes) const;

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
    
    // Initialize Gauss points
    void initializeGaussPoints();

private:
    // Gauss point data
    std::vector<double> m_gaussR, m_gaussS, m_gaussT, m_gaussW;
    
    // Cached Jacobian data
    mutable std::vector<Matrix3d> m_jacobianInverse;
    mutable std::vector<double> m_jacobianDet;
};

} // namespace RgFem

#endif // RGHEX20ELEMENT_H