#pragma once
#include "FileImport.h"
#include <FECore/FESurfacePairConstraint.h>

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEFileSection
{
protected:
	//! missing primary surface
	class MissingPrimarySurface : public FEFileException
	{
	public:
		MissingPrimarySurface();
	};

	//! missing secondary surface
	class MissingSecondarySurface : public FEFileException
	{
	public:
		MissingSecondarySurface();
	};

public:
	FEBioContactSection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseLinearConstraint     (XMLTag& tag);

protected:
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};

class FEBioContactSection4 : public FEBioContactSection
{
public:
	FEBioContactSection4(FEFileImport* im) : FEBioContactSection(im) {}
	void Parse(XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci);
};
