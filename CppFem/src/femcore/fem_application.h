/*****************************************************************//**
 * \file   fem_application.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef FEM_APP_HH
#define FEM_APP_HH

#include "problemmode.h"
#include <memory>

namespace fem
{
	class EngngModel;

	class FemApp
	{
	public:
		std::unique_ptr<EngngModel> InstanciateProblem(DataReader& dr, ProblemMode mode, int contextFlag, EngngModel* _master, bool parallelFlag);
	};
} //end namespace fem

#endif
