#include "FERestartImport.h"
#include "femcore/Solver/FESolver.h"
#include "femcore/FEAnalysis/FEAnalysis.h"
#include "femcore/Solver/FESolver.h"
#include "femcore/FEModel.h"
#include "basicio/DumpFile.h"
#include "femcore/TimeStep/FETimeStepController.h"
#include "FEBioLoadDataSection.h"
#include "FEBioStepSection.h"
//#include "FEBioStepSection3.h"
#include "FEBioStepSection4.h"

void FERestartControlSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	++tag;
	do
	{
		if      (tag == "time_steps"        ) tag.value(pstep->m_ntime);
		else if (tag == "final_time"        ) tag.value(pstep->m_final_time);
		else if (tag == "step_size"         ) tag.value(pstep->m_dt0);
		else if (tag == "time_stepper")
		{
			if (pstep->m_timeController == nullptr) pstep->m_timeController = RANGO_NEW<FETimeStepController>( & fem,"");
			FETimeStepController& tc = *pstep->m_timeController;
			++tag;
			do
			{
				if      (tag == "max_retries") tag.value(tc.m_maxretries);
				else if (tag == "opt_iter"   ) tag.value(tc.m_iteopt);
				else if (tag == "dtmin"      ) tag.value(tc.m_dtmin);
				else throw XMLReader::InvalidTag(tag);

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "plot_level")
		{
			char szval[256];
			tag.value(szval);
			if      (strcmp(szval, "PLOT_NEVER"        ) == 0) pstep->SetPlotLevel(FE_PLOT_NEVER);
			else if (strcmp(szval, "PLOT_MAJOR_ITRS"   ) == 0) pstep->SetPlotLevel(FE_PLOT_MAJOR_ITRS);
			else if (strcmp(szval, "PLOT_MINOR_ITRS"   ) == 0) pstep->SetPlotLevel(FE_PLOT_MINOR_ITRS);
			else if (strcmp(szval, "PLOT_MUST_POINTS"  ) == 0) pstep->SetPlotLevel(FE_PLOT_MUST_POINTS);
			else if (strcmp(szval, "PLOT_FINAL"        ) == 0) pstep->SetPlotLevel(FE_PLOT_FINAL);
			else if (strcmp(szval, "PLOT_STEP_FINAL"   ) == 0) pstep->SetPlotLevel(FE_PLOT_STEP_FINAL);
			else if (strcmp(szval, "PLOT_AUGMENTATIONS") == 0) pstep->SetPlotLevel(FE_PLOT_AUGMENTATIONS);
			else throw XMLReader::InvalidValue(tag);
		}
		else if (tag == "plot_stride")
		{
			int n = 1;
			tag.value(n);
			if (n < 1) throw XMLReader::InvalidValue(tag);
			pstep->SetPlotStride(n);
		}
		else if (tag == "solver")
		{
			FEAnalysis* step = fem.GetCurrentStep();
			FESolver* solver = step->GetFESolver();
			if (solver == nullptr) throw XMLReader::InvalidTag(tag);

			++tag;
			do
			{
				if (ReadParameter(tag, solver) == false)
				{
					throw XMLReader::InvalidTag(tag);
				}
				++tag;
			} while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	// we need to reevaluate the time step size and end time
	pstep->m_dt = pstep->m_dt0;
//	pstep->m_tend = pstep->m_tstart = pstep->m_ntime*pstep->m_dt0;

}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FERestartImport::FERestartImport()
{
	m_newSteps = 0;
}

FERestartImport::~FERestartImport()
{
	
}

int FERestartImport::StepsAdded() const
{
	return m_newSteps;
}

//-----------------------------------------------------------------------------
bool FERestartImport::Load(FEModel& fem, const char* szfile)
{
	// open the XML file
	if (m_xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening restart file %s\n", szfile);

	m_builder = new FEModelBuilder(fem);

	m_szdmp[0] = 0;
	m_newSteps = 0;

	m_map["Control" ] = new FERestartControlSection(this);

	// loop over child tags
	bool ret = true;
	try
	{
		// find the root element
		XMLTag tag;
		if (m_xml.FindTag("febio_restart", tag) == false) return errf("FATAL ERROR: File does not contain restart data.\n");

		// check the version number
		const char* szversion = tag.AttributeValue("version");
		int nversion = -1;
		if      (strcmp(szversion, "1.0") == 0) nversion = 1;
		else if (strcmp(szversion, "2.0") == 0) nversion = 2;
		else if (strcmp(szversion, "3.0") == 0) nversion = 3;
		else if (strcmp(szversion, "4.0") == 0) nversion = 4;

		if (nversion == -1) return errf("FATAL ERROR: Incorrect restart file version\n");

		// Add the Step section for version 2
		if (nversion == 2)
		{
			// set the file version to make sure we are using the correct format
			SetFileVerion(0x0205);

			// make sure we can redefine curves in the LoadData section
			FEBioLoadDataSection* lcSection = new FEBioLoadDataSection(this);
			lcSection->SetRedefineCurvesFlag(true);
			m_map["LoadData"] = lcSection;

			m_map["Step"] = new FEBioStepSection25(this);
		}

		// Add the Step section for version 3
		if (nversion == 3)
		{
			// set the file version to make sure we are using the correct format
			SetFileVerion(0x0300);

			// make sure we can redefine curves in the LoadData section
			FEBioLoadDataSection3* lcSection = new FEBioLoadDataSection3(this);
			lcSection->SetRedefineCurvesFlag(true);
			m_map["LoadData"] = lcSection;

			//m_map["Step"] = new FEBioStepSection3(this);
		}

		// Add the Step section for version 4
		if (nversion == 4)
		{
			// set the file version to make sure we are using the correct format
			SetFileVerion(0x0400);

			// make sure we can redefine curves in the LoadData section
			FEBioLoadDataSection3* lcSection = new FEBioLoadDataSection3(this);
			lcSection->SetRedefineCurvesFlag(true);
			m_map["LoadData"] = lcSection;

			m_map["Step"] = new FEBioStepSection4(this);
		}

		// the first section has to be the archive
		++tag;
		if (tag != "Archive") return errf("FATAL ERROR: The first element must be the archive name\n");
		char szar[256];
		tag.value(szar);

		// open the archive
		DumpFile ar(fem);
		if (ar.Open(szar) == false) return errf("FATAL ERROR: failed opening restart archive\n");

		// read the archive
		fem.Serialize(ar);

		// set the module name
		GetBuilder()->SetActiveModule(fem.GetModuleName());

		// keep track of the number of steps
		int steps0 = fem.Steps();

		// read the rest of the restart input file
		ret = ParseFile(tag);

		// count nr of newly added steps
		int steps1 = fem.Steps();
		m_newSteps = steps1 - steps0;
	}
	catch (XMLReader::Error& e)
	{
		fprintf(stderr, "FATAL ERROR: %s\n", e.what());
		ret = false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: unrecoverable error (line %d)\n", m_xml.GetCurrentLine());
		ret = false;
	}

	// close the XML file
	m_xml.Close();

	// we're done!
	return ret;
}
