#include "stdafx.h"
#include "FEBioInitialSection3.h"
#include <femcore/FEInitialCondition.h>

FEBioInitialSection3::FEBioInitialSection3(FEFileImport* pim) : FEFileSection(pim) 
{
}

void FEBioInitialSection3::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if (tag == "ic") ParseIC(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

void FEBioInitialSection3::ParseIC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();

	// read the type attribute
	const char* sztype = tag.AttributeValue("type");

	// try to allocate the initial condition
	FEInitialCondition* pic = fecore_new<FEInitialCondition>(sztype, fem);
	if (pic == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// add it to the model
	GetBuilder()->AddInitialCondition(pic);

	FENodalIC* nic = dynamic_cast<FENodalIC*>(pic);
	if (nic)
	{
		// read required node_set attribute
		const char* szset = tag.AttributeValue("node_set");
		FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
		if (nodeSet == nullptr) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);
		nic->SetNodeSet(nodeSet);
	}

	// Read the parameter list
	ReadParameterList(tag, pic);
}
