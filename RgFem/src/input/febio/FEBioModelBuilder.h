#pragma once
#include "input/FEBio/FEModelBuilder.h"

class FEBioModel;

class FEBioModelBuilder : public FEModelBuilder
{
public:
	FEBioModelBuilder(FEModel& fem);

public:
	FEDomain* CreateDomain(FE_Element_Spec espec, FEMaterial* mat) override;
	void AddMaterial(FEMaterial* mat) override;
	void AddRigidComponent(FEStepComponent* prc) override;
};

