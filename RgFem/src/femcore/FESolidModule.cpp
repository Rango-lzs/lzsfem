#include "femcore/FESolidModule.h"
#include "femcore/DOFS.h"
#include "femcore/FEModel.h"

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

static const char* GetVariableName(MECH_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT:
            return "displacement";
            break;
        case SHELL_ROTATION:
            return "shell rotation";
            break;
        case RIGID_ROTATION:
            return "rigid rotation";
            break;
        case SHELL_DISPLACEMENT:
            return "shell displacement";
            break;
        case VELOCTIY:
            return "velocity";
            break;
        case SHELL_VELOCITY:
            return "shell velocity";
            break;
        case SHELL_ACCELERATION:
            return "shell acceleration";
            break;
    }

    assert(false);
    return nullptr;
}

FESolidModule::FESolidModule() {}

void FESolidModule::InitModel(FEModel* fem)
{
	// Allocate degrees of freedom
	DOFS& dofs = fem->GetDOFS();
	int varD = dofs.AddVariable(GetVariableName(DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable(GetVariableName(SHELL_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable(GetVariableName(RIGID_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varV = dofs.AddVariable(GetVariableName(VELOCTIY), VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
	int varSU = dofs.AddVariable(GetVariableName(SHELL_DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varSU, 0, "sx");
	dofs.SetDOFName(varSU, 1, "sy");
	dofs.SetDOFName(varSU, 2, "sz");
	int varSV = dofs.AddVariable(GetVariableName(SHELL_VELOCITY), VAR_VEC3);
	dofs.SetDOFName(varSV, 0, "svx");
	dofs.SetDOFName(varSV, 1, "svy");
	dofs.SetDOFName(varSV, 2, "svz");
	int varSA = dofs.AddVariable(GetVariableName(SHELL_ACCELERATION), VAR_VEC3);
	dofs.SetDOFName(varSA, 0, "sax");
	dofs.SetDOFName(varSA, 1, "say");
	dofs.SetDOFName(varSA, 2, "saz");
}
