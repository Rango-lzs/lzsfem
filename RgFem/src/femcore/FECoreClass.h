#pragma once
#include "FEObjectBase.h"
#include "FEM_EXPORT.h"

// Generic base class for classes that don't fit in the super-class structure
class FEM_EXPORT FECoreClass : public FEObjectBase
{
	FECORE_SUPER_CLASS(FECLASS_ID)

public:
	FECoreClass(FEModel* fem = nullptr);

	DECLARE_FECORE_CLASS();
};
