#ifndef RGHEX8ELEMENT_H
#define RGHEX8ELEMENT_H

#include "RgLinearSolid3dElement.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include <array>
#include <vector>
#include <string>

// forward declarations for project types
class RgMaterial;
class Matrix3ds;

namespace RgFem {

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
class RgHex8Element : public RgLinearSolid3dElement {
public:
    static constexpr int kNodeCount = 8;

    // Constructors and Destructors
    RgHex8Element();
    explicit RgHex8Element(bool fullInt);
    RgHex8Element(const RgHex8Element& other);
    RgHex8Element& operator=(const RgHex8Element& other);
    virtual ~RgHex8Element();

    //使用ElementTraitsStore来管理ElementTraits,
    static RgElementTraits* fullIntTraits();
    static RgElementTraits* reduceIntTraits();

    // Stiffness and Mass Matrices
    virtual void calculateStiffnessMatrix(Matrix& K) const override;
    virtual void calculateMassMatrix(Matrix& M) const override;
    virtual void calculateDampingMatrix(Matrix& C) const override;
     
    // Internal force vector (residual)
    virtual void calculateInternalForceVector(std::vector<double>& F) const override;

    // Strain and Stress Calculations
    virtual void calculateStress(FEMaterialPoint& matPt, StressTensor& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, StrainTensor& strain) override;
    
    // Strain energy calculation
    virtual double calculateStrainEnergy() const override;
    // Kinetic energy calculation
    virtual double calculateKineticEnergy() const override;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;
   
protected:
    // Compute B matrix (strain-displacement matrix) at natural coordinates
    void computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const;
};

} // namespace RgFem

#endif // RGHEX8ELEMENT_H