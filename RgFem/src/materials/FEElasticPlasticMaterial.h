#pragma once
#include "materials/FEElasticMaterial.h"
#include "FEElasticPlasticMaterialPoint.h"
#include "femcore/Domain/FEDomainParameter.h"

//-----------------------------------------------------------------------------
//! Base class for elastic-plastic materials

class FEM_EXPORT FEElasticPlasticMaterial : public FEElasticMaterial
{
    DECLARE_META_CLASS(FEElasticPlasticMaterial, FEElasticMaterial);

public:
	//! constructor 
	FEElasticPlasticMaterial();

	//! destructor
	~FEElasticPlasticMaterial();

	//! create material point data for this material
	FEMaterialPointData* CreateMaterialPointData() override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! calculate Cauchy stress
    Matrix3ds Stress(FEMaterialPoint& pt) override;

    //! calculate spatial tangent stiffness
    tens4ds Tangent(FEMaterialPoint& pt) override;

    //! calculate the 2nd P-K stress
    Matrix3ds PK2Stress(FEMaterialPoint& pt, const Matrix3ds E) override;

public:
    // Material parameters
    FEParamDouble m_E;      //!< Young's modulus
    FEParamDouble m_v;      //!< Poisson's ratio
    FEParamDouble m_sy;     //!< Yield stress

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FEElasticPlasticStress : public FEDomainParameter
{
public:
	FEElasticPlasticStress();
	FEParamValue value(FEMaterialPoint& mp) override;
};