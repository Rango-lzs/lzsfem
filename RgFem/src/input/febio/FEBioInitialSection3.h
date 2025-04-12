#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection3 : public FEFileSection
{
public:
	FEBioInitialSection3(FEFileImport* pim);
	void Parse(XMLTag& tag);
	void ParseIC(XMLTag& tag);
};
