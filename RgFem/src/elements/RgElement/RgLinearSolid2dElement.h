#pragma once

#include "RgSolid2dElement.h"

using RgMatrix = Matrix;
using RgVector = Vector;

//! RgLinearSolid2dElement - Linear 2D solid element
//! Used for small strain, small displacement analysis (geometrically linear)
class FEM_EXPORT RgLinearSolid2dElement : public RgSolid2dElement
{
public:
    //! default constructor
    RgLinearSolid2dElement() = default;

    //! copy constructor
    RgLinearSolid2dElement(const RgLinearSolid2dElement& el) = default;

    //! assignment operator
    RgLinearSolid2dElement& operator=(const RgLinearSolid2dElement& el) = default;

    //! destructor
    virtual ~RgLinearSolid2dElement() = default;


    // Stiffness and Mass Matrices
    virtual void calculateStiffnessMatrix(RgMatrix& K) override;
    virtual void calculateMassMatrix(RgMatrix& M) override;
    virtual void calculateDampingMatrix(RgMatrix& C) override;

    // Internal force vector (residual)
    virtual void calculateInternalForceVector(std::vector<double>& F) override;

    // Strain and Stress Calculations
    virtual void getStress(RgMaterialPoint& matPt, StressTensor& stress);
    virtual void getStrain(RgMaterialPoint& matPt, StrainTensor& strain);

    void updateStress(RgMaterialPoint& matPt);
    void updateStrain(RgMaterialPoint& matPt);

    // Energy calculations
    virtual double calculateStrainEnergy() override;
    virtual double calculateKineticEnergy() override;
};