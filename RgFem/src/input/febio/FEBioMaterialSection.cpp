#include "FEBioMaterialSection.h"
#include "femcore/FEModel.h"
#include "materials/FEMaterial.h"
#include "logger/log.h"
#include <sstream>

//-----------------------------------------------------------------------------
//! This function creates a material by checking the type attribute against
//! registered materials. Also, if the tag defines attributes (other than
//! type and name), the material is offered a chance to process the attributes.
FEMaterial* FEBioMaterialSection::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);
	
	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// create a new material of this type
	FEMaterial* pmat = GetBuilder()->CreateMaterial(sztype);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	return pmat;
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		if (tag == "material")
		{
			// check that the ID attribute is defined and that it 
			// equals the number of materials + 1.
			int nid = -1;
			tag.AttributeValue("id", nid);
			int nmat = fem.Materials();
			if (nid != nmat+1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// make sure that the name is unique
			std::string name;
			const char* szname = tag.AttributeValue("name", true);
			if (szname == nullptr)
			{
				std::stringstream ss;
				ss << "Material" << nid;
				name = ss.str();

				feLogWarningEx((&fem), "Material %d has no name.\nIt was given the name %s.", nid, name.c_str());
			}
			else name = szname;

			FEMaterial* mat = fem.FindMaterial(name);
			if (mat)
			{
				throw XMLReader::InvalidAttributeValue(tag, "name");
			}

			// create a material from this tag
			FEMaterial* pmat = CreateMaterial(tag); assert(pmat);
			pmat->SetName(name);

			// set the material's ID
			++m_nmat;
			pmat->SetID(m_nmat);

			// parse the material parameters
			ReadParameterList(tag, pmat);

			// add the material
			GetBuilder()->AddMaterial(pmat);
		}
		else throw XMLReader::InvalidTag(tag);

		// read next tag
		++tag;
	}
	while (!tag.isend());
}

//===============================================================================

//-----------------------------------------------------------------------------
//! This function creates a material by checking the type attribute against
//! registered materials. Also, if the tag defines attributes (other than
//! type and name), the material is offered a chance to process the attributes.
FEMaterial* FEBioMaterialSection3::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);

	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// create a new material of this type
	FEMaterial* pmat = GetBuilder()->CreateMaterial(sztype);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	return pmat;
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
void FEBioMaterialSection3::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		if (tag == "material")
		{
			// check that the ID attribute is defined and that it 
			// equals the number of materials + 1.
			int nid = -1;
			tag.AttributeValue("id", nid);
			int nmat = fem.Materials();
			if (nid != nmat + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// make sure that the name is unique
			const char* szname = tag.AttributeValue("name");
			FEMaterial* mat = fem.FindMaterial(szname);
			if (mat)
			{
				throw XMLReader::InvalidAttributeValue(tag, "name");
			}

			// create a material from this tag
			FEMaterial* pmat = CreateMaterial(tag); assert(pmat);
			pmat->SetName(szname);

			// add the material
			fem.AddMaterial(pmat);
			++m_nmat;

			// set the material's ID
			pmat->SetID(m_nmat);

			// parse the material parameters
			ReadParameterList(tag, pmat);
		}
		else throw XMLReader::InvalidTag(tag);

		// read next tag
		++tag;
	} while (!tag.isend());
}
