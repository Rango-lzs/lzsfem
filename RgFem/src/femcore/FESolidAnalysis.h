#pragma once
#include "femcore/FEAnalysis/FEAnalysis.h"
class FEM_EXPORT FESolidAnalysis : public FEAnalysis
{
	DECLARE_META_CLASS(FESolidAnalysis, FEAnalysis);
public:
	enum SolidAnalysisType {
		STATIC,
		DYNAMIC
	};

public:
	FESolidAnalysis();

	DECLARE_PARAM_LIST();
};
