#include "RgDomainMap.h"
#include "FEModel.h"
#include "RgSolidDomain.h"
#include "RgShellDomain.h"
#include "RgBeamDomain.h"
#include "RgTrussDomain.h"
#include "FESolidDomainFactory.h"

//-----------------------------------------------------------------------------
RgDomainMap::RgDomainMap(FEModel* pfem) : RgDomainList(pfem)
{
}

//-----------------------------------------------------------------------------
RgDomainMap::~RgDomainMap()
{
}

//-----------------------------------------------------------------------------
RgDomain* RgDomainMap::CreateDomain(const char* sztype)
{
	FEModel& fem = *GetFEModel();
	RgDomain* pdom = 0;
/*
	if      (strcmp(sztype, "solid") == 0) pdom = new FESolidDomain(&fem);
	else if (strcmp(sztype, "shell") == 0) pdom = new FEShellDomain(&fem);
	else if (strcmp(sztype, "beam" ) == 0) pdom = new FEBeamDomain (&fem);
	else if (strcmp(sztype, "truss") == 0) pdom = new FETrussDomain(&fem);
	else if (strncmp(sztype, "solid-", 6) == 0)
	{
		// see if we can create it using the factory
		pdom = FESolidDomainFactory::CreateSolidDomain(sztype+6, &fem);
	}
*/
	return pdom;
}

//-----------------------------------------------------------------------------
void RgDomainMap::AddDomain(RgDomain* pd)
{
	Add(pd);
}