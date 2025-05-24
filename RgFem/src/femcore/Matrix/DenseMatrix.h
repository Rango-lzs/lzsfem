#pragma once
#include "SparseMatrix.h"


//=============================================================================
//! This class implements a full Matrix
//! that is a Matrix that stores all its elements.

class FEM_EXPORT DenseMatrix : public SparseMatrix
{
public:
	// con/de-structor
	DenseMatrix();
	~DenseMatrix();

	// create a Matrix of particular size
	void Create(int rows, int cols);

	// retrieve Matrix data
	double& operator () (int i, int j) { return m_pr[i][j]; }

public:
	// zero Matrix elements
	void Zero() override;

	// create a Matrix from a spares Matrix profile
	void Create(SparseMatrixProfile& mp) override;

	// assemble Matrix into sparse Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lm) override;

	//! assemble a Matrix into the sparse Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	// clear all data
	void Clear() override;

	bool check(int i, int j) override { return true; }
	void add(int i, int j, double v) override { m_pr[i][j] += v; }
	void set(int i, int j, double v) override { m_pr[i][j] = v; }
	double diag(int i) override { return m_pr[i][i]; }

protected:
	double*		m_pd;	//!< Matrix values
	double**	m_pr;	//!< pointers to rows
};

