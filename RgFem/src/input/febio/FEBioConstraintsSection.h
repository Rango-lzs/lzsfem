#pragma once
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Constraints Section (format 2.5)
class FEBioConstraintsSection25 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection25(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
};
