#pragma once
#include "FEBioImport.h"

class FEBioLoadsSection3 : public FEFileSection
{
public:
	FEBioLoadsSection3(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	//void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
};
