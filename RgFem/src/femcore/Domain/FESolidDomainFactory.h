#pragma once

#include "femcore/FEDomainFatory.h"
//-----------------------------------------------------------------------------
class FESolidDomainFactory : public FEDomainFactory
{
public:
	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
