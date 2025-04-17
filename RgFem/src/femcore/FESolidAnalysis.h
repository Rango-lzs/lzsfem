#pragma once
#include "femcore/FEAnalysis/FEAnalysis.h"
class FEM_EXPORT FESolidAnalysis : public FEAnalysis
{
public:
	enum SolidAnalysisType {
		STATIC,
		DYNAMIC
	};

public:
	FESolidAnalysis(FEModel* fem);

	DECLARE_PARAM_LIST();
};
