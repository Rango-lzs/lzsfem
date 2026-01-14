#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEDomainMap;
class FESurfaceMap;
class FENodeDataMap;

//-----------------------------------------------------------------------------
// MeshData Section (introduced in febio_spec 2.5)
class FEBioMeshDataSection : public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[RgElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseShellThickness(XMLTag& tag, RgElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, RgElementSet& set);
	void ParseMaterialAxes  (XMLTag& tag, RgElementSet& set);
	void ParseMaterialData  (XMLTag& tag, RgElementSet& set, const std::string& name);
	void ParseMaterialFiberProperty(XMLTag& tag, RgElementSet& set);

private:
	void ParseElementData(XMLTag& tag, RgElementSet& set, std::vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseDataArray(XMLTag& tag, FEDataArray& map, const char* sztag);
	RgElement* ParseElement(XMLTag& tag, RgElementSet& set);

private:
	std::vector<RgElement*> m_pelem;
};

//-----------------------------------------------------------------------------
// MeshData Section for febio_spec 3.0
class FEBioMeshDataSection3 : public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[RgElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection3(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalData(XMLTag& tag);
	void ParseSurfaceData(XMLTag& tag);
	void ParseElementData(XMLTag& tag);

protected:
//	void ParseModelParameter(XMLTag& tag, FEParamValue param);
//	void ParseRgMaterialPointData(XMLTag& tag, FEParamValue param);

protected:
	void ParseShellThickness(XMLTag& tag, RgElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, RgElementSet& set);
	void ParseMaterialAxes(XMLTag& tag, RgElementSet& set);
	void ParseMaterialAxesProperty(XMLTag& tag, RgElementSet& set);

private:
	void ParseElementData(XMLTag& tag, RgElementSet& set, std::vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseSurfaceData(XMLTag& tag, FESurfaceMap& map);
	void ParseNodeData   (XMLTag& tag, FENodeDataMap& map);
};

//-----------------------------------------------------------------------------
// MeshData Section for febio_spec 4.0
class FEBioMeshDataSection4: public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[RgElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection4(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalData(XMLTag& tag);
	void ParseSurfaceData(XMLTag& tag);
	void ParseElementData(XMLTag& tag);

private:
	void ParseNodeData(XMLTag& tag, FENodeDataMap& map);
	void ParseSurfaceData(XMLTag& tag, FESurfaceMap& map);
	void ParseElementData(XMLTag& tag, RgElementSet& set, std::vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseShellThickness(XMLTag& tag, RgElementSet& set);
	void ParseMaterialAxes(XMLTag& tag, RgElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, RgElementSet& set);
};
