#pragma once
#include "materials/FEMaterial.h"
#include "datastructure/tens4d.h"

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FEM_EXPORT FESolidMaterial : public FEMaterial
{
public:
	//! constructor
	FESolidMaterial(FEModel* pfem);

	//柯西应力 causy stress
	virtual Matrix3ds Stress(FEMaterialPoint& pt) = 0;

	//calculate spatial tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//calculate the 2nd Piola-Kirchhoff stress at material point
	virtual Matrix3ds PK2Stress(FEMaterialPoint& pt, const Matrix3ds E);

	//calculate material tangent stiffness at material poin
	virtual tens4dmm MaterialTangent(FEMaterialPoint& pt, const Matrix3ds E);

    // calculate secant tangent stiffness at material point 割线刚度
    virtual tens4dmm SecantTangent(FEMaterialPoint& pt, bool mat = false);

	//the material density
	void SetDensity(double d);

	//evaluate density
	virtual double Density(FEMaterialPoint& pt);

	//! Is this a rigid material or not
	virtual bool IsRigid() const { return false; }

	//返回切线或割线刚度阵
	tens4dmm SolidTangent(FEMaterialPoint& pt);

	virtual Matrix3ds SecantStress(FEMaterialPoint& pt, bool PK2 = false);
	virtual bool UseSecantTangent() { return false; }

protected:
	FEParamDouble m_density;	//!< material density

	DECLARE_PARAM_LIST();
};
