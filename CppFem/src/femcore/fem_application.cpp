/*****************************************************************//**
 * \file   fem_application.cpp
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#include "fem_application.h"
#include "engngm.h"
#include "fem_factory.h"
#include "inputrecord.h"
#include "datareader.h"
#include "error.h"

#include <cstring>
#if defined ( __GNUC__ ) && defined ( HAVE_EXECINFO_H )
#include <cxxabi.h>
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>
#endif

namespace fem 
{
	std::unique_ptr<EngngModel> FemApp::InstanciateProblem(DataReader& dr, ProblemMode mode, int contextFlag, EngngModel* _master, bool parallelFlag)
	{
		std::string problemName, dataOutputFileName, desc;

		dataOutputFileName = dr.giveOutputFileName();
		desc = dr.giveDescription();

		/* here we need copy of input record. The pointer returned by dr.giveInputRecord can (and will)
		 * be updated as reading e-model components (nodes, etc). But we need this record being available
		 * through the whole e-model instanciation
		 */
		auto emodelir = dr.giveInputRecord(DataReader::IR_emodelRec, 1).clone();
		emodelir->giveRecordKeywordField(problemName); ///@todo Make this function robust, it can't be allowed to fail (the record keyword is not a normal field-id)

		auto problem = FemFactory::instance().createEngngModel(problemName.c_str(), 1, _master);
		if (!problem) {
			FEM_WARNING("Failed to construct engineering model of type \"%s\".\n", problemName.c_str());
			return NULL;
		}

		problem->setProblemMode(mode);
		problem->setParallelMode(parallelFlag);

		if (contextFlag) {
			problem->setContextOutputMode(COM_Always);
		}

		problem->instanciateYourself(dr, *emodelir, dataOutputFileName.c_str(), desc.c_str());
		emodelir->finish();

		return problem;
	}
} // end namespace fem
