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

	//! calculate stress at material point
	virtual Matrix3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate the 2nd Piola-Kirchhoff stress at material point
	virtual Matrix3ds PK2Stress(FEMaterialPoint& pt, const Matrix3ds E);

	//! calculate material tangent stiffness at material point
	virtual tens4dmm MaterialTangent(FEMaterialPoint& pt, const Matrix3ds E);

    //! calculate secant tangent stiffness at material point
    virtual tens4dmm SecantTangent(FEMaterialPoint& pt, bool mat = false);

	//! return the material density
	void SetDensity(const double d);

	//! evaluate density
	virtual double Density(FEMaterialPoint& pt);

	//! Is this a rigid material or not
	virtual bool IsRigid() const { return false; }

	tens4dmm SolidTangent(FEMaterialPoint& pt);

	virtual Matrix3ds SecantStress(FEMaterialPoint& pt, bool PK2 = false);
	virtual bool UseSecantTangent() { return false; }

protected:
	FEParamDouble m_density;	//!< material density
};
