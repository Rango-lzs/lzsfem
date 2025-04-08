#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Material Section
class FEM_EXPORT FEBioMaterialSection : public FEFileSection
{
public:
	FEBioMaterialSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};

//-----------------------------------------------------------------------------
// Material Section
class FEM_EXPORT FEBioMaterialSection3 : public FEFileSection
{
public:
	FEBioMaterialSection3(FEFileImport* pim) : FEFileSection(pim) { m_nmat = 0; }
	void Parse(XMLTag& tag);

protected:
	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};
