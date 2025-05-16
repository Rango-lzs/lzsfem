#pragma once
#include "femcore/fem_export.h"
//-----------------------------------------------------------------------------
//! Creation of domains are a little more elaborate and deviate from the usual
//! factory methods.
class FEM_EXPORT FEDomainFactory
{
public:
    FEDomainFactory()
    {
    }
    virtual ~FEDomainFactory()
    {
    }

    virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat) = 0;
};

//-----------------------------------------------------------------------------
class FESolidDomainFactory : public FEDomainFactory
{
public:
	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
