#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Step Section (4.0 format)
class FEBioStepSection4 : public FEFileSection
{
public:
	FEBioStepSection4(FEFileImport* pim);
	void Parse(XMLTag& tag);
};
