#pragma once
#include "input/febio/FEBioImport.h"
#include "femcore/FERgModel.h"

//-----------------------------------------------------------------------------
// Mesh section
class FEBioMeshSection4 : public FEBioFileSection
{
public:
	FEBioMeshSection4(FEBioImport* pim);

	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag, FEBModel::Part* part);
	void ParseSurfaceSection    (XMLTag& tag, FEBModel::Part* part);
	void ParseElementSection    (XMLTag& tag, FEBModel::Part* part);
	void ParseNodeSetSection    (XMLTag& tag, FEBModel::Part* part);
	void ParseElementSetSection (XMLTag& tag, FEBModel::Part* part);
	void ParsePartListSection   (XMLTag& tag, FEBModel::Part* part);
	void ParseEdgeSection       (XMLTag& tag, FEBModel::Part* part);
	void ParseSurfacePairSection(XMLTag& tag, FEBModel::Part* part);
	void ParseDiscreteSetSection(XMLTag& tag, FEBModel::Part* part);
};
