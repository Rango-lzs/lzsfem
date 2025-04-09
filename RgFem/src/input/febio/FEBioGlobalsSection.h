#pragma once
#include "FileImport.h"

//-----------------------------------------------------------------------------
//! Globals Section
class FEM_EXPORT FEBioGlobalsSection : public FEFileSection
{
public:
	FEBioGlobalsSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag) override;

protected:
	void ParseConstants   (XMLTag& tag);
	void ParseGlobalData  (XMLTag& tag);
	void ParseVariables   (XMLTag& tag);
};
