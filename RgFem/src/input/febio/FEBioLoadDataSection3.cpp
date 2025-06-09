#include "input/febio/FEBioLoadDataSection.h"
#include "femcore/FEModel.h"
#include "femcore/FELoadController.h"
#include "femcore/FEPointFunction.h"

//-----------------------------------------------------------------------------
FEBioLoadDataSection3::FEBioLoadDataSection3(FEFileImport* pim) : FEFileSection(pim) 
{
	m_redefineCurves = false;
}

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection3::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	
	++tag;
	do
	{
		if (tag == "load_controller")
		{
			// load curve ID
			int nid;
			tag.AttributeValue("id", nid);

			// get the type
			const char* sztype = tag.AttributeValue("type");

			// get the number of load curves
			int nlc = fem.LoadControllers();

			// create the controller
            FELoadController* plc = RANGO_NEW<FELoadController>(&fem, sztype);
			if (plc == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type");

			// set the ID
			plc->SetID(nid - 1);

			// see if this refers to a valid curve
			if (m_redefineCurves && ((nid > 0) && (nid <= nlc)))
			{
				fem.ReplaceLoadController(nid - 1, plc);
			}
			else
			{
				// check that the ID is one more than the number of load curves defined
				// This is to make sure that the ID's are in numerical order and no values are skipped.
				if (nid != nlc + 1)
				{
					delete plc;
					throw XMLReader::InvalidAttributeValue(tag, "id");
				}

				// add the controller
				fem.AddLoadController(plc);
			}

			// read the parameter list
			ReadParameterList(tag, plc);

			if (m_redefineCurves)
			{
				// We only get here during restart, in which case the new load controllers
				// will not get a chance to be initialized, so we'll do it here.
				plc->Init();
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
