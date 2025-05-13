#pragma once
#include "FEBioImport.h"
#include "FEBModel.h"

//-----------------------------------------------------------------------------
// MeshDomains section
class FEBioMeshDomainsSection4 : public FEBioFileSection
{
public:
	FEBioMeshDomainsSection4(FEBioImport* pim);

	void Parse(XMLTag& tag);

protected:
	void ParseSolidDomainSection(XMLTag& tag);
	void ParseShellDomainSection(XMLTag& tag);
	void ParseBeamDomainSection(XMLTag& tag);

private:
	void BuildNLT();

private:
	std::vector<int>	m_NLT;
	int					m_noff;
};
