#ifndef RGBEAM3DELEMENT_H
#define RGBEAM3DELEMENT_H

#include "RgBeamElement.h"
#include <array>
#include <vector>

namespace RgFem {

// forward declarations
class Matrix3d;
class Vector3d;
class Matrix3ds;

/*
    RgBeam3dElement
    - 2-node linear Timoshenko beam element in 3D
    - 1D isoparametric element with 6 DOFs per node (3 translations + 3 rotations)
    - Derives from RgBeamElement
    - Features:
      * Linear shape functions for displacement
      * Accounts for shear deformation (Timoshenko formulation)
      * 2-point Gauss quadrature integration
      * Supports small strain/small rotation analysis
      * Node distribution: 2 nodes at element ends
      * Total DOF: 12 (6 per node: ux, uy, uz, rx, ry, rz)
      * Full 3D bending in any orientation
*/
class RgBeam3dElement : public RgBeamElement
{
public:
    static constexpr int kNodeCount = 2;
    static constexpr int kDofsPerNode = 6;  // ux, uy, uz, rx, ry, rz
    static constexpr int kTotalDofs = kNodeCount * kDofsPerNode;
    static constexpr int kNumGaussPoints = 2;

    // Constructors and Destructors
    RgBeam3dElement();
    explicit RgBeam3dElement(const std::array<int, kNodeCount>& nodeIds);
    RgBeam3dElement(const RgBeam3dElement& other);
    RgBeam3dElement& operator=(const RgBeam3dElement& other);
    virtual ~RgBeam3dElement();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementClass() const;

    // Element Properties
    virtual int getNumberOfNodes() const override { return kNodeCount; }
    virtual int getNumberOfGaussPoints() const { return kNumGaussPoints; }
    virtual int getDofsPerNode() const { return kDofsPerNode; }
    virtual int getTotalDofs() const { return kTotalDofs; }

    // Initialize element traits
    virtual void initTraits() override;

    // Geometric operations
    virtual Vector3d evaluateCoordinates(double naturalCoord) const;
    virtual double evaluateLength() const;
    Vector3d getAxis() const;

    // Shape function evaluations
    void evaluateShapeFunctions(double r, std::vector<double>& N) const;
    void evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const;

    // FEM Matrix Calculations
    virtual void calculateStiffnessMatrix(Matrix& K) const override;
    virtual void calculateMassMatrix(Matrix& M) const override;
    virtual void calculateDampingMatrix(Matrix& C) const override;

    // Strain and Stress Calculations
    virtual void calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain) override;

    // Loading Operations
    virtual void applyBodyForce(const Vector3d& force, Vector& F) const;
    virtual void applyDistributedLoad(const Vector3d& load, Vector& F) const;
    virtual void applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const;
    virtual void applyMoment(int nodeId, const Vector3d& moment, Vector& F) const;

    // Material point access
    virtual FEMaterialPoint* getMaterialPoint(int gaussPtId);
    virtual const FEMaterialPoint* getMaterialPoint(int gaussPtId) const;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

    // Utility functions
    double getCharacteristicLength() const;
    double getVolume() const;

    // Beam specific properties
    double getElementLength() const { return evaluateLength(); }
    void setMomentOfInertia(double Iy, double Iz) { m_Iy = Iy; m_Iz = Iz; }
    void setArea(double A) { m_A = A; }
    void setShearFactor(double alpha) { m_shearFactor = alpha; }
    void setTorsionalConstant(double Ip) { m_Ip = Ip; }
    void setOrientationVector(const Vector3d& v) { m_orient = v; }

    // Geometric nonlinear analysis (large rotation)
    void updateRotations(const double* displacements);
    Matrix3d computeRotationMatrix() const;

protected:
    // Compute transformation matrix (local to global)
    Matrix3d computeTransformationMatrix() const;

    // Compute element stiffness in local coordinates
    void computeLocalStiffness(Matrix& K_local) const;

    // Initialize Gauss points
    void initializeGaussPoints();

private:
    // Gauss point data
    std::vector<double> m_gaussR, m_gaussW;

    // Beam properties
    double m_Iy;  // Second moment of inertia about y-axis
    double m_Iz;  // Second moment of inertia about z-axis
    double m_Ip;  // Polar moment of inertia (torsional constant)
    double m_A;   // Cross-sectional area
    double m_shearFactor;  // Shear correction factor (default 5/6 for rectangular)

    // Element orientation
    Vector3d m_orient;  // Orientation vector (if not along x-axis)
    
    // For large rotation analysis
    std::vector<Matrix3d> m_rotation;  // Rotation matrix at each Gauss point
};

} // namespace RgFem

#endif // RGBEAM3DELEMENT_H
