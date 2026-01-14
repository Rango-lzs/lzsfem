#pragma once
#include "materials/FESolidMaterial.h"
#include "FEElasticMaterialPoint.h"
#include "femcore/Domain/FEDomainParameter.h"

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEM_EXPORT FEElasticMaterial : public FESolidMaterial
{
    DECLARE_META_CLASS(FEElasticMaterial, FESolidMaterial);

public:
	//! constructor 
	FEElasticMaterial();

	//! destructor
	~FEElasticMaterial();

	//! create material point data for this material
	RgMaterialPointData* CreateRgMaterialPointData() override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(RgMaterialPoint& pt);
    
    // get the elastic material
    virtual FEElasticMaterial* GetElasticMaterial() { return this; }

public:
	//! evaluates approximation to Cauchy stress using forward difference
	Matrix3ds SecantStress(RgMaterialPoint& pt, bool PK2 = false) override;

public:
    virtual double StrongBondSED(RgMaterialPoint& pt) { return StrainEnergyDensity(pt); }
    virtual double WeakBondSED(RgMaterialPoint& pt) { return 0; }


	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FEElasticStress : public FEDomainParameter
{
public:
	FEElasticStress();
	FEParamValue value(RgMaterialPoint& mp) override;
};
