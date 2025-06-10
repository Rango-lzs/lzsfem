#pragma once
#include "femcore/FEModule.h"

enum MECH_VARIABLE
{
    DISPLACEMENT,
    SHELL_ROTATION,
    RIGID_ROTATION,
    SHELL_DISPLACEMENT,
    VELOCTIY,
    SHELL_VELOCITY,
    SHELL_ACCELERATION,
};

const char* GetVariableName(MECH_VARIABLE var);

class FEM_EXPORT FESolidModule : public FEModule
{
public:
	FESolidModule();
	void InitModel(FEModel* fem) override;
};
