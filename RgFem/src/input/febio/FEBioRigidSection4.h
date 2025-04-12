#pragma once
#include "FEBioImport.h"

class FEBioRigidSection4 : public FEFileSection
{
public:
	FEBioRigidSection4(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidBC(XMLTag& tag);
	void ParseRigidIC(XMLTag& tag);
	void ParseRigidLoad(XMLTag& tag);
	void ParseRigidConnector(XMLTag& tag);
};
