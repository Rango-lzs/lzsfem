#pragma once
#include "elements/RgElement/RgElement.h"
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix3d.h"
#include <vector>
#include "elements/RgGaussPoint.h"

//! Forward declarations
class RgSolidElementTraits;

namespace RgFem
{
    class NaturalCoord;
}

//! Corresponds to Abaqus Continuum(Solid) elements, which are multidimensional (3D)
class FEM_EXPORT RgSolidElement : public RgElement
{
public:
	//! default constructor
	RgSolidElement() = default;
	~RgSolidElement() = default;

	RgSolidElement(const RgSolidElement& el);
	RgSolidElement& operator = (const RgSolidElement& el);

	//! return spatial dimension (3 for 3D solids)
	virtual int dim() = 0;
    virtual int dofs() const;

    // n : the n-th gauss point
    virtual RgFem::RgGaussPoint gaussPoint(int n) const;

    virtual std::vector<double> evalH(int n);
    // [3,N] (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(int n);
    //[6,N] (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(int n);

    //! 这些接口计算任意点的形函数
    virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord);
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord);
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord);
 
	virtual void calculateStiffnessMatrix(Matrix& K) const override;         // Km = ∫ B^T D B dV (default linear assembly)
	virtual void calculateMassMatrix(Matrix& M) const override;
	virtual void calculateDampingMatrix(Matrix& C) const override;
	virtual void calculateInternalForceVector(std::vector<double>& F) const override;

	virtual void computeBMatrix(const Vector3d& xi, Matrix& B) const = 0;

    virtual void calculateStress(FEMaterialPoint& matPt, StressTensor& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, StrainTensor& strain) override;

    virtual double calculateStrainEnergy() const;
    virtual double calculateKineticEnergy() const;
    virtual void Serialize(DumpStream& ar) override;

    virtual Matrix3d evaluateJacobian(const Vector3d& naturalCoord) const;
    double evaluateJacobianDeterminant(const Vector3d& naturalCoord) const;
    virtual Matrix3d evaluateJacobianInverse(const Vector3d& naturalCoord) const;

protected:
	virtual void computeBMatrixAtGauss(int gp, Matrix& B) const;
	virtual void getConstitutiveTangentAtGauss(int gp, Matrix& D) const;
};