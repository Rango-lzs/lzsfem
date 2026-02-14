#pragma once

#include "RgStructureElement.h"
#include <vector>

// forward declarations
class Matrix3d;
class Vector3d;

/*
    RgBeamElement
    - Abstract base class for beam elements
    - Derived from RgStructureElement
    - Provides common interface for 1D structural elements
    - Features:
      * 1D isoparametric element framework
      * Small strain/small rotation analysis
      * Support for body forces and concentrated loads
      * Material point integration
*/
class FEM_EXPORT RgBeamElement : public RgStructureElement
{
public:
    // Constructors and Destructors
    RgBeamElement();
    RgBeamElement(const RgBeamElement& other);
    RgBeamElement& operator=(const RgBeamElement& other);
    virtual ~RgBeamElement();

    // Element Type Identification
    virtual ElementCategory elementCategory() const override { return ElementCategory::FE_ELEM_BEAM; }

    // Element Properties
    virtual int getNumberOfGaussPoints() const = 0;
    virtual int getDofsPerNode() const = 0;
    virtual int getTotalDofs() const = 0;

    // Geometric operations
    virtual double evaluateLength() const = 0;
    virtual Vector3d getAxis() const = 0;

    // Shape function evaluations
    virtual void evaluateShapeFunctions(double r, std::vector<double>& N) const = 0;
    virtual void evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const = 0;

    // FEM Matrix Calculations
    virtual void calculateStiffnessMatrix(Matrix& K) override;
    virtual void calculateMassMatrix(Matrix& M) override ;
    virtual void calculateDampingMatrix(Matrix& C) override;

    //// Strain and Stress Calculations
    //virtual void calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress) override;
    //virtual void calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain) override;

    // Loading Operations
    virtual void applyBodyForce(const Vector3d& force, Vector& F) const = 0;
    virtual void applyDistributedLoad(const Vector3d& load, Vector& F) const = 0;
    virtual void applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const = 0;

    // Material point access
    virtual FEMaterialPoint* getMaterialPoint(int gaussPtId) = 0;
    virtual const FEMaterialPoint* getMaterialPoint(int gaussPtId) const = 0;

    // Utility functions
    virtual double getCharacteristicLength() const { return evaluateLength(); }
    virtual double getVolume() const = 0;

    // Serialization
    virtual void Serialize(DumpStream& ar) override = 0;

protected:
    // Initialize Gauss points (called by derived classes)
    virtual void initializeGaussPoints() = 0;
};

