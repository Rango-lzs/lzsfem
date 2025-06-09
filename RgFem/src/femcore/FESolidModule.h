#pragma once
#include "femcore/FEModule.h"

class FEM_EXPORT FESolidModule : public FEModule
{
public:
	FESolidModule();
	void InitModel(FEModel* fem) override;
};
