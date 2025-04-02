#pragma once

#include "matrix.h"
#include "vector.h"
#include "LinearSolver.h"
#include "FENewtonStrategy.h"

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class FEM_EXPORT BFGSSolver : public FENewtonStrategy
{
public:
	//! constructor
	BFGSSolver(FEModel* fem);

	//! New initialization method
	bool Init() override;

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) override;

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b) override;

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver
	int				m_neq;		//!< number of equations

	// BFGS update vectors
	matrix			m_V;		//!< BFGS update vector
	matrix			m_W;		//!< BFGS update vector
	vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors

	vector<double>	tmp;

	DECLARE_FECORE_CLASS();
};
