#pragma once
#include "FEBioBoundarySection.h"

//-----------------------------------------------------------------------------
// version 3.0
class FEBioBoundarySection3 : public FEBioBoundarySection
{
public:
	FEBioBoundarySection3(FEFileImport* imp) : FEBioBoundarySection(imp) {}
	void Parse(XMLTag& tag);
	void ParseBC(XMLTag& tag);
	void ParseBCRigid(XMLTag& tag);	// read rigid contact section
	void ParseLinearConstraint(XMLTag& tag);
};
