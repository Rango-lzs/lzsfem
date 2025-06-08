#pragma once

#include "datastructure/Matrix.h"
#include "LinearSolver.h"
#include "FENewtonStrategy.h"

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class FEM_EXPORT BFGSSolver : public FENewtonStrategy
{
    DECLARE_META_CLASS(BFGSSolver, FENewtonStrategy);

public:
	//! constructor
	BFGSSolver();

	//! New initialization method
	bool Init() override;

	//! perform a BFGS udpate
	bool Update(double s, std::vector<double>& ui, std::vector<double>& R0, std::vector<double>& R1) override;

	//! solve the equations
	void SolveEquations(std::vector<double>& x, std::vector<double>& b) override;

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver
	int				m_neq;		//!< number of equations

	// BFGS update vectors
	Matrix			m_V;		//!< BFGS update std::vector
	Matrix			m_W;		//!< BFGS update std::vector
	std::vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors

	std::vector<double>	tmp;

	DECLARE_PARAM_LIST();
};
