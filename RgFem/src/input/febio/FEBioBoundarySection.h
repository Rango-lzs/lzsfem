#pragma once
#include "FEBioImport.h"
#include "femcore/FESurfacePairConstraint.h"
#include <map>

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEFileSection
{
public:
	FEBioBoundarySection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseBCFix         (XMLTag& tag);
	void ParseBCPrescribe   (XMLTag& tag);
	void ParseContactSection(XMLTag& tag);
	void ParseConstraints   (XMLTag& tag);
	void ParseSpringSection (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* psi);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);

protected:
	void BuildNodeSetMap();

	void AddFixedBC(FENodeSet* set, int bc);

protected:
	std::map<std::string, FENodeSet*>	m_NodeSet;	// map for faster lookup of node sets
};

