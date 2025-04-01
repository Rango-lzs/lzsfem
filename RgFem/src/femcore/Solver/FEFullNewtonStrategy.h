#pragma once

#include "LinearSolver.h"
#include "FENewtonStrategy.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEM_EXPORT FEFullNewtonStrategy : public FENewtonStrategy
{
public:
	//! constructor
	FEFullNewtonStrategy(FEModel* fem);

	//! New initialization method
	bool Init() override;

	//! perform a BFGS udpate
    bool Update(double s, std::vector<double>& ui, std::vector<double>& R0, std::vector<double>& R1) override;

	//! solve the equations
    void SolveEquations(std::vector<double>& x, std::vector<double>& b) override;

public:
	// keep a pointer to the linear solver
	LinearSolver* m_plinsolve;	//!< pointer to linear solver
};
