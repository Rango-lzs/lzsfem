#pragma once
#include "LinearSolver.h"
#include "femcore/Matrix/DenseMatrix.h"

//-----------------------------------------------------------------------------
//! LU decomposition solver

//! This solver performs an LU decomposition and uses a backsolving algorithm
//! to solve the equations.
//! This solver uses the FullMatrix class and therefore is not the preferred
//! solver. It should only be used for small problems and only when the other
//! solvers are not adequate.

class FEM_EXPORT LUSolver : public LinearSolver
{
public:
	//! constructor
	LUSolver(FEModel* fem = nullptr);

	//! Pre-process data
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! solve using factored matrix
	bool BackSolve(double* x, double* b) override;

	//! Clean-up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(MatrixType ntype) override;

	//! Set the matrix
	void SetMatrix(FECore::DenseMatrix* pA);

protected:
	std::vector<int>		indx;	//!< indices
	FECore::DenseMatrix*	m_pA;	//!< sparse matrix
};
