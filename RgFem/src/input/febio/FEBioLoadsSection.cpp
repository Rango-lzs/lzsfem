#include "FEBioLoadsSection.h"
#include "femcore/FEBodyLoad.h"
#include "femcore/FEModel.h"
#include "femcore/FENodalLoad.h"
#include "femcore/FESurfaceLoad.h"
#include "femcore/FEScalarValuator.h"
//#include "femcore/FEEdgeLoad.h"
//#include "femcore/FEEdge.h"

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "nodal_load"  ) ParseNodalLoad  (tag);
		else if (tag == "surface_load") ParseSurfaceLoad(tag);
		//else if (tag == "edge_load"   ) ParseEdgeLoad   (tag);
		else if (tag == "body_load"   ) ParseBodyLoad   (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
	const char* sztype = tag.AttributeValue("type");
    FEBodyLoad* pbl = RANGO_NEW<FEBodyLoad>(&fem, sztype);
	if (pbl == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// see if a name was defined
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pbl->SetName(szname);

	// see if a specific domain was referenced
	const char* szpart = tag.AttributeValue("elem_set", true);
	if (szpart)
	{
		FEMesh& mesh = fem.GetMesh();
		RgElementSet* elset = mesh.FindElementSet(szpart);
		if (elset == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szpart);
		pbl->SetDomainList(elset);
	}

	ReadParameterList(tag, pbl);
	fem.AddModelLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseNodalLoad(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type
	const char* sztype = tag.AttributeValue("type");

	// create nodal load
    FENodalLoad* pfc = RANGO_NEW<FENodalLoad>(GetFEModel() , sztype);
	if (pfc == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type");

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pfc->SetName(szname);

	// get the node set
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// assign the node set
	pfc->SetNodeSet(nodeSet);

	// add it to the model
	GetBuilder()->AddNodalLoad(pfc);

	// read parameters
	ReadParameterList(tag, pfc);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create surface load
	std::string typeString = tag.AttributeValue("type");
    FESurfaceLoad* psl = RANGO_NEW<FESurfaceLoad>(&fem, typeString.c_str());
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// read required surface attribute
	const char* surfaceName = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(surfaceName);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", surfaceName);

	// create a surface from this facet set
    FESurface* psurf = RANGO_NEW<FESurface>(&fem, "");
	GetBuilder()->BuildSurface(*psurf, *pface);

	// assign it
	mesh.AddSurface(psurf);
	psl->SetSurface(psurf);

	// add it to the model
	GetBuilder()->AddSurfaceLoad(psl);

	// read the parameters
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (ReadParameter(tag, psl) == false)
			{
				// some special handling for fluidnormalvelocity
				if ((tag == "value") && (strcmp(typeString.c_str(), "fluid normal velocity") == 0))
				{
					// The value parameter no longer exists, but for backward compatibility
					// we map it to the "velocity" parameter.

					// get the surface data attribute
					const char* szsurfdata = tag.AttributeValue("surface_data");

					// find the surface map
					FEDataMap* surfMap = mesh.FindDataMap(szsurfdata);
					if (surfMap == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface_data");

					// get the velocity parameter
					FEParam* p = psl->GetParameter("velocity");
					FEParamDouble& v = p->value<FEParamDouble>();

					// we'll map the velocity value to the map's scale factor
					double s = 1.0;
					if (v.isConst()) s = v.constValue(); else throw XMLReader::InvalidTag(tag);

					// create map and assign
					FEMappedValue* map = new FEMappedValue();
					map->setDataMap(surfMap);
					map->setScaleFactor(s);
					v.setValuator(map);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			++tag;
		} 
		while (!tag.isend());
	}
}

////-----------------------------------------------------------------------------
//void FEBioLoadsSection3::ParseEdgeLoad(XMLTag& tag)
//{
//	FEModel& fem = *GetFEModel();
//	FEMesh& mesh = fem.GetMesh();
//
//	// create edge load
//	const char* sztype = tag.AttributeValue("type");
//	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(sztype, &fem);
//	if (pel == 0) throw XMLReader::InvalidTag(tag);
//
//	// create a new edge
//	FEEdge* pedge = new FEEdge(&fem);
//	mesh.AddEdge(pedge);
//	pel->SetEdge(pedge);
//
//	// get the segment set
//	const char* szedge = tag.AttributeValue("edge");
//	FESegmentSet* pset = mesh.FindSegmentSet(szedge);
//	if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "edge", szedge);
//	if (GetBuilder()->BuildEdge(*pedge, *pset) == false) throw XMLReader::InvalidTag(tag);
//
//	// add edge load to model
//	GetBuilder()->AddEdgeLoad(pel);
//
//	// read the parameters
//	ReadParameterList(tag, pel);
//}
