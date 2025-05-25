#pragma once
#include "FileImport.h"
#include "input/XML/XMLReader.h"

//-----------------------------------------------------------------------------
class FERestartControlSection : public FEFileSection
{
public:
	FERestartControlSection(FEFileImport* reader) : FEFileSection(reader) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! Restart input file reader.
class FEBIOXML_API FERestartImport : public FEFileImport
{
public:
	FERestartImport();
	virtual ~FERestartImport();

	bool Load(FEModel& fem, const char* szfile);

	int StepsAdded() const;

public:
	char		m_szdmp[256];	// user defined restart file name

protected:
	XMLReader	m_xml;			// the file reader
	int		m_newSteps;		// nr of new steps added
};
