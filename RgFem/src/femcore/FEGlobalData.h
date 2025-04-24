#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! This class can be used to define global model data and will be placed in the
//! global date section of the FEModel class
class FEM_EXPORT FEGlobalData : public FEModelComponent
{
    DECLARE_META_CLASS(FEGlobalData, FEModelComponent);

public:
	//! constructor
	FEGlobalData(FEModel* fem);
};
