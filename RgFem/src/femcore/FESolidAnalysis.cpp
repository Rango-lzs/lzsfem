#include "femcore/FESolidAnalysis.h"

DEFINE_META_CLASS(FESolidAnalysis, FEAnalysis, "solid");

BEGIN_PARAM_DEFINE(FESolidAnalysis, FEAnalysis)
	// The analysis parameter is already defined in the FEAnalysis base class. 
	// Here, we just need to set the enum values for the analysis parameter.
	FindParameterFromData(&m_nanalysis)->setEnums("STATIC\0DYNAMIC\0");
END_PARAM_DEFINE()

FESolidAnalysis::FESolidAnalysis() : FEAnalysis()
{

}
