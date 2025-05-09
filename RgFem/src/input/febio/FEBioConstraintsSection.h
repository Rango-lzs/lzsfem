#pragma once
#include "FEBioImport.h"

class FEFacetSet;

// (Base class. Don't use this directly!)
class FEBioConstraintsSection : public FEFileSection
{
public:
    FEBioConstraintsSection(FEFileImport* pim)
        : FEFileSection(pim)
    {
    }

protected:
    bool ParseSurfaceSection(XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};
//-----------------------------------------------------------------------------
// Constraints Section (format 2.5)
class FEBioConstraintsSection25 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection25(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
};
