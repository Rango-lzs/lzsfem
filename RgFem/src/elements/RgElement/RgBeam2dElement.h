#ifndef RGBEAM2DELEMENT_H
#define RGBEAM2DELEMENT_H

#include "RgLinearBeamElement.h"
#include <array>
#include <vector>



// forward declarations
class Matrix3d;
class Vector3d;

/*
    RgBeam2dElement
    - 2-node linear Timoshenko beam element
    - 1D isoparametric element with 3 DOFs per node (2 displacements + 1 rotation)
    - Derives from RgLinearBeamElement
    - Features:
      * Linear shape functions for displacement
      * Accounts for shear deformation (Timoshenko formulation)
      * 2-point Gauss quadrature integration
      * Supports small strain/small rotation analysis
      * Node distribution: 2 nodes at element ends
      * Total DOF: 6 (3 per node)
*/
class RgBeam2dElement : public RgLinearBeamElement
{
public:
    static constexpr int kNodeCount = 2;
    static constexpr int kDofsPerNode = 3;  // ux, uy, rz
    static constexpr int kTotalDofs = kNodeCount * kDofsPerNode;
    static constexpr int kNumGaussPoints = 2;

    // Constructors and Destructors
    RgBeam2dElement();
    explicit RgBeam2dElement(const std::array<int, kNodeCount>& nodeIds);
    RgBeam2dElement(const RgBeam2dElement& other);
    RgBeam2dElement& operator=(const RgBeam2dElement& other);
    virtual ~RgBeam2dElement();

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
    void setMomentOfInertia(double I) { m_I = I; }
    void setArea(double A) { m_A = A; }
    void setShearFactor(double alpha) { m_shearFactor = alpha; }

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
    double m_I;  // Second moment of inertia
    double m_A;  // Cross-sectional area
    double m_shearFactor;  // Shear correction factor (default 5/6 for rectangular)
};



#endif // RGBEAM2DELEMENT_H
