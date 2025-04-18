#pragma once
#include "femcore/FEGlobalVector.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEMechModel;


//The FEResidualVector implements a global vector that stores the residual.
class FEM_EXPORT FEResidualVector : public FEGlobalVector
{
public:
	FEResidualVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr);
	~FEResidualVector();

	//Assemble the element vector into this global vector
	void Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom = false) override;
	//Assemble into this global vector
	void Assemble(int node, int dof, double f) override;
};
