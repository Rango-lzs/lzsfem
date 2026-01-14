#pragma once

#include "RgSolid3dElement.h"

    //! RgLinearSolid3dElement - Linear 3D solid element
    //! Used for small strain, small displacement analysis (geometrically linear)
    //! Features:
    //!   - Linear strain-displacement relationship
    //!   - Standard B-matrix formulation
    //!   - Material stiffness integration via Gauss quadrature
    //!   - Compatible with isotropic and anisotropic materials
class FEM_EXPORT RgLinearSolid3dElement : public RgSolid3dElement
{
public:
    //! default constructor
    RgLinearSolid3dElement() = default;

    //! copy constructor
    RgLinearSolid3dElement(const RgLinearSolid3dElement& el) = default;

    //! assignment operator
    RgLinearSolid3dElement& operator=(const RgLinearSolid3dElement& el) = default;

    //! destructor
    virtual ~RgLinearSolid3dElement() = default;

    // Stiffness and Mass Matrices
    virtual void calculateStiffnessMatrix(Matrix& K) override;
    virtual void calculateMassMatrix(Matrix& M) override;
    virtual void calculateDampingMatrix(Matrix& C) override;

    // Internal force vector (residual)
    virtual void calculateInternalForceVector(std::vector<double>& F) override;

    // Strain and Stress Calculations
    virtual void getStress(RgMaterialPoint& matPt, StressTensor& stress);
    virtual void getStrain(RgMaterialPoint& matPt, StrainTensor& strain);

    void updateStress(RgMaterialPoint& matPt);
    void updateStrain(RgMaterialPoint& matPt);

    // Strain energy calculation
    virtual double calculateStrainEnergy() override;
    // Kinetic energy calculation
    virtual double calculateKineticEnergy() override;
};