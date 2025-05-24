#pragma once

#include "LinearSolver.h"
#include "femcore/Matrix/SkylineMatrix.h"

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class FEM_EXPORT SkylineSolver : public LinearSolver
{
public:
	//! constructor
	SkylineSolver(FEModel* fem);

	//! Preprocess 
	bool PreProcess() override;

	//! Factor Matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(double* x, double* b) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse Matrix
	SparseMatrix* CreateSparseMatrix(MatrixType ntype) override;

private:
	SkylineMatrix*	m_pA;
};
