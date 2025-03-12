#pragma once
#include "FEModelComponent.h"
#include "FEMaterialPoint.h"
#include "FEModelParam.h"
#include "FEDomainList.h"
#include "FEDomainParameter.h"

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEDomain;
class DumpStream;

//-----------------------------------------------------------------------------
class FEM_EXPORT FEMaterialBase : public FEModelComponent
{
public:
	FEMaterialBase(FEModel* fem);

	//! returns a pointer to a new material point object
	virtual FEMaterialPointData* CreateMaterialPointData();

	//! Update specialized material points at each iteration
	virtual void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp);

	// evaluate local coordinate system at material point
	virtual mat3d GetLocalCS(const FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types
//! From this class all other material classes are derived.

class FEM_EXPORT FEMaterial : public FEMaterialBase
{
	FECORE_SUPER_CLASS(FEMATERIAL_ID)
	FECORE_BASE_CLASS(FEMaterial)

public:
	FEMaterial(FEModel* fem);
	virtual ~FEMaterial();

	//! performs initialization
	bool Init() override;

	//! get a domain parameter
	FEDomainParameter* FindDomainParameter(const std::string& paramName);

	// evaluate local coordinate system at material point
	mat3d GetLocalCS(const FEMaterialPoint& mp) override;

	// set the (local) material axis valuator
	void SetMaterialAxis(FEMat3dValuator* val);

protected:
	FEMat3dValuator*	m_Q;			//!< local material coordinate system

public:
	//! Assign a domain to this material
	void AddDomain(FEDomain* dom);

	//! get the domaint list
	FEDomainList& GetDomainList() { return m_domList; }

protected:
	void AddDomainParameter(FEDomainParameter* p);

private:
	FEDomainList	m_domList;		//!< list of domains that use this material

	std::vector<FEDomainParameter*>	m_param;	//!< list of domain variables

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Material properties are classes that can only be defined as properties of other materials
class FEM_EXPORT FEMaterialProperty : public FEMaterialBase
{
	FECORE_SUPER_CLASS(FEMATERIALPROP_ID)

public:
	FEMaterialProperty(FEModel* fem);

	// evaluate local coordinate system at material point
	mat3d GetLocalCS(const FEMaterialPoint& mp) override;

	DECLARE_FECORE_CLASS();
};
