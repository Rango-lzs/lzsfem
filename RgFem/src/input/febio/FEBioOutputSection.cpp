#include "FEBioOutputSection.h"
#include "basicio/NodeDataRecord.h"
#include "basicio/FaceDataRecord.h"
#include "basicio/ElementDataRecord.h"
#include "basicio/ObjectDataRecord.h"
#include "basicio/NLConstraintDataRecord.h"
#include "basicio/SurfaceDataRecord.h"
#include "basicio/DomainDataRecord.h"
#include "basicio/FEModelDataRecord.h"
#include "femcore/FEModel.h"
#include "femcore/FSPath.h"
#include "femcore/FEPlotDataStore.h"
#include "femcore/FESurface.h"
#include "basicio/DataRecord.h"
//#include "femcore/FESurfaceDataRecord.h"

bool string_to_int_vector(const char* szlist, std::vector<int>& list)
{
	list.clear();
	if (szlist == nullptr) return false;
	if (szlist[0] == 0) return true;

	char* ch = nullptr;
	char* sz = (char*)szlist;
	do
	{
		ch = strchr(sz, ',');
		if (ch) *ch = 0;
		int n0, n1, nn;
		int nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 0: return false;
			break;
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
		}

		for (int i = n0; i <= n1; i += nn) list.push_back(i);

		if (ch) *ch = ',';
		sz = ch + 1;
	}
	while (ch != 0);

	return true;
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "logfile" ) ParseLogfile(tag);
		else if (tag == "plotfile") ParsePlotfile(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParseLogfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// Get the feb file path
	const char* szpath = GetFileReader()->GetFilePath();

	// see if the log file has any attributes
	const char* szlog = tag.AttributeValue("file", true);

	// If the log filename is not a path, but just a filename, make the filename
	// relative to the input file, otherwise, leave the path alone.
	if(szlog)
	{
		if(!FSPath::isPath(szlog))
		{
			char szfile[1024] = {0};
			sprintf(szfile, "%s%s", szpath, szlog);
			GetFEBioImport()->SetLogfileName(szfile);
		}
		else
		{
			GetFEBioImport()->SetLogfileName(szlog);
		}

	}

	if (tag.isleaf()) return;
	++tag;
	do
	{
		// get the (optional) file attribute
		char szfilename[1024] = { 0 };
		const char* szfile = tag.AttributeValue("file", true);
		if (szfile)
		{
			// if we have a path, prepend the path's name
			if (szpath && szpath[0])
			{
				sprintf(szfilename, "%s%s", szpath, szfile);
			}
			else strcpy(szfilename, szfile);
			szfile = szfilename;
		}

		// get other attributes
		const char* szdelim = tag.AttributeValue("delim", true);
		const char* szformat = tag.AttributeValue("format", true);

		bool bcomment = true;
		const char* szcomment = tag.AttributeValue("comments", true);
		if (szcomment != 0)
		{
			if (strcmp(szcomment, "on") == 0) bcomment = true;
			else if (strcmp(szcomment, "off") == 0) bcomment = false;
		}

		// get the data attribute
		const char* szdata = tag.AttributeValue("data");

		// get the name attribute
		const char* szname = tag.AttributeValue("name", true);

		// allocate data record
		DataRecord* pdr = nullptr;
		if (tag == "node_data")
		{
            pdr = RANGO_NEW<DataRecord>(&fem,"node_data");

			const char* sztmp = "set";
			if (GetFileReader()->GetFileVersion() >= 0x0205) sztmp = "node_set";
			const char* sz = tag.AttributeValue(sztmp, true);
			if (sz)
			{
				FENodeSet* pns = mesh.FindNodeSet(sz);
				if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, sztmp, sz);

				std::vector<int> items;
				int n = pns->Size();
				assert(n);
				items.resize(n);
				for (int i = 0; i < n; ++i) items[i] = (*pns)[i] + 1;

				pdr->SetItemList(items);
			}
			else
			{
				std::vector<int> items;
				string_to_int_vector(tag.szvalue(), items);
				pdr->SetItemList(items);
			}
		}
		else if (tag == "face_data")
		{
            pdr = RANGO_NEW<DataRecord>(&fem ,"face_data");

			const char* sz = tag.AttributeValue("surface");
			FESurface* surf = mesh.FindSurface(sz);
			if (surf == nullptr)
			{
				FEFacetSet* pfs = mesh.FindFacetSet(sz);
				if (pfs == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);

				surf = new FESurface(&fem);
				surf->Create(*pfs);
				surf->SetName(sz);
				surf->Init();
				mesh.AddSurface(surf);
			}

			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);

			// TODO: This is a bit of a hack, because the face data record needs an FEItemList, but FESurface is derived from that.
			FEFacetSet* fset = surf->GetFacetSet();
			fset->SetSurface(surf);
			pdr->SetItemList(fset, items);
		}
		else if (tag == "element_data")
		{
            pdr = RANGO_NEW<DataRecord>(&fem ,"element_data");

			const char* sztmp = "elset";
			if (GetFileReader()->GetFileVersion() >= 0x0205) sztmp = "elem_set";

			const char* sz = tag.AttributeValue(sztmp, true);
			if (sz)
			{
				std::vector<int> dummy;
				FEElementSet* pes = mesh.FindElementSet(sz);
				if (pes == 0) throw XMLReader::InvalidAttributeValue(tag, sztmp, sz);
				pdr->SetItemList(pes, dummy);
			}
			else
			{
				std::vector<int> items;
				string_to_int_vector(tag.szvalue(), items);
				pdr->SetItemList(items);
			}
		}
		else if (tag == "rigid_body_data")
		{
            pdr = RANGO_NEW<DataRecord>(&fem ,"rigid_body_data");
			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);
			pdr->SetItemList(items);
		}
        else if (tag == "rigid_connector_data")
        {
            pdr = RANGO_NEW<DataRecord>(&fem ,"rigid_connector_data");
			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);
			pdr->SetItemList(items);
        }
        else if (tag == "surface_data")
        {
            FESurfaceDataRecord* prec = new FESurfaceDataRecord(&fem);
			pdr = prec;
			const char* sz = tag.AttributeValue("surface");
			if (sz)
			{
				int surfIndex = mesh.FindSurfaceIndex(sz);
				if (surfIndex == -1) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
				prec->SetSurface(surfIndex);
			}
        }
        else if (tag == "domain_data")
        {
            FEDomainDataRecord* prec = new FEDomainDataRecord(&fem);
			pdr = prec;
			const char* sz = tag.AttributeValue("domain");
			if (sz)
			{
				int domainIndex = mesh.FindDomainIndex(sz);
				if (domainIndex == -1) throw XMLReader::InvalidAttributeValue(tag, "domain", sz);
				prec->SetDomain(domainIndex);
			}
        }
        else if (tag == "model_data")
        {
            pdr = new FEModelDataRecord(&fem);
        }
		else throw XMLReader::InvalidTag(tag);

		if (pdr)
		{
			pdr->SetData(szdata);
			if (szname != 0) pdr->SetName(szname); else pdr->SetName(szdata);
			if (szfile) pdr->SetFileName(szfile);
			if (szdelim != 0) pdr->SetDelim(szdelim);
			if (szformat != 0) pdr->SetFormat(szformat);
			pdr->SetComments(bcomment);
			GetFEBioImport()->AddDataRecord(pdr);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParsePlotfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	FEPlotDataStore& plotData = fem.GetPlotDataStore();

	// get the plot file type. Must be "febio"!
	const char* sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if ((strcmp(sz, "febio") != 0) && (strcmp(sz, "febio2") != 0)) throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	else sz = "febio";
	plotData.SetPlotFileType(sz);

	// get the optional plot file name
	const char* szplt = tag.AttributeValue("file", true);

	// If the plot filename is not a path, but just a filename, make the filename
	// relative to the input file, otherwise, leave the path alone.
	if(szplt)
	{
		if(!FSPath::isPath(szplt))
		{
			// Get the feb file path
			const char* szpath = GetFileReader()->GetFilePath();

			char szfile[1024] = {0};
			sprintf(szfile, "%s%s", szpath, szplt);
			GetFEBioImport()->SetPlotfileName(szfile);
		}
		else
		{
			GetFEBioImport()->SetPlotfileName(szplt);
		}

	}

	// read and store the plot variables
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (tag == "var")
			{
				// get the variable name
				const char* szt = tag.AttributeValue("type");

				// get the item list
				std::vector<int> item;
				if (tag.isempty() == false) tag.value(item);

                // see if a surface is referenced
                const char* szsurf = tag.AttributeValue("surface", true);
                const char* szeset = tag.AttributeValue("elem_set", true);
                if (szsurf)
                {
                    // make sure this tag does not have any children
                    if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);
                    
                    // see if we can find the facet set
                    FEMesh& m = GetFEModel()->GetMesh();
					FEFacetSet* ps = m.FindFacetSet(szsurf);
                    
                    // create a surface from the facet set
                    if (ps)
                    {
                        // create a new surface
                        FESurface* psurf = RANGO_NEW < FESurface>(&fem, "");
                        fem.GetMesh().AddSurface(psurf);
                        if (GetBuilder()->BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);

                        // Add the plot variable
                        const std::string& surfName = psurf->GetName();
						plotData.AddPlotVariable(szt, item, surfName.c_str());
                    }
                    else throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);
                }
				else if (szeset)
				{
					// make sure this tag does not have any children
					if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

					// see if we can find the facet set
					FEMesh& m = GetFEModel()->GetMesh();
					FEElementSet* ps = m.FindElementSet(szeset);

					// create a surface from the facet set
					if (ps)
					{
						// Add the plot variable
						plotData.AddPlotVariable(szt, item, szeset);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "elem_set", szeset);
				}
                else
                {
                    // Add the plot variable
					plotData.AddPlotVariable(szt, item);
                }
			}
			else if (tag=="compression")
			{
				int ncomp;
				tag.value(ncomp);
				plotData.SetPlotCompression(ncomp);
			}
			++tag;
		}
		while (!tag.isend());
	}
}
