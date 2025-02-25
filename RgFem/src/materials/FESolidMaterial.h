#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/tens4d.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FEBIOMECH_API FESolidMaterial : public FEMaterial
{
public:
	//! constructor
	FESolidMaterial(FEModel* pfem);

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate the 2nd Piola-Kirchhoff stress at material point
	virtual mat3ds PK2Stress(FEMaterialPoint& pt, const mat3ds E);

	//! calculate material tangent stiffness at material point
	virtual tens4dmm MaterialTangent(FEMaterialPoint& pt, const mat3ds E);

    //! calculate secant tangent stiffness at material point
    virtual tens4dmm SecantTangent(FEMaterialPoint& pt, bool mat = false);

	//! return the material density
	void SetDensity(const double d);

	//! evaluate density
	virtual double Density(FEMaterialPoint& pt);

	//! Is this a rigid material or not
	virtual bool IsRigid() const { return false; }

	tens4dmm SolidTangent(FEMaterialPoint& pt);

	virtual mat3ds SecantStress(FEMaterialPoint& pt, bool PK2 = false);
	virtual bool UseSecantTangent() { return false; }

protected:
	FEParamDouble	m_density;	//!< material density

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FESolidMaterial)
};
