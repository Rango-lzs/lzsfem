#pragma once
#include "femcore/FEPrescribedDOF.h"

class FEPrescribedDisplacement : public FEPrescribedDOF
{
    DECLARE_META_CLASS(FEPrescribedDisplacement, FEPrescribedDOF);

public:
	FEPrescribedDisplacement();

private:
	DECLARE_PARAM_LIST();
};
