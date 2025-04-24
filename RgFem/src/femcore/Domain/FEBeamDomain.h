#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for beam domains
class FEM_EXPORT FEBeamDomain : public FEDomain
{
    DECLARE_META_CLASS(FEBeamDomain, FEDomain);

public:
	FEBeamDomain(FEModel* pm);
};
