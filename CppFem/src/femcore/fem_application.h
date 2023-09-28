/*****************************************************************//**
 * \file   fem_application.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef FEM_APP_HH
#define FEM_APP_HH

#include <memory>

namespace fem
{
	//”√¥¶≤ªœÍ
	enum ProblemMode { 
		_processor,
		_postProcessor
	};

	/// Corresponds to macro- and micro-problem in multiscale simulations.
	enum ProblemScale {
		macroScale,
		microScale
	};

	class EngngModel;

	class FemApp
	{
	public:
		std::unique_ptr<EngngModel> InstanciateProblem(DataReader& dr, ProblemMode mode, int contextFlag, EngngModel* _master, bool parallelFlag);
	};
} //end namespace fem

#endif
