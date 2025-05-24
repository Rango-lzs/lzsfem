#pragma once
#include "SparseMatrix.h"

//=============================================================================
//! Implements a sparse Matrix using the skyline storage

//! This class implements a symmetric sparse Matrix where only the values
//! below the skyline are stored.

class FEM_EXPORT SkylineMatrix : public SparseMatrix
{
public:
	SkylineMatrix();
	virtual ~SkylineMatrix();

public: // from SparseMatrix

	void Zero() override;

	void Clear() override;

	void Create(SparseMatrixProfile& mp) override;

	void Assemble(const Matrix& ke, const std::vector<int>& lm) override;

	//! assemble a Matrix into the sparse Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	void add(int i, int j, double v) override;

	void set(int i, int j, double v) override;

	// NOTE: This is not implemented yet!
	bool check(int i, int j) override;

	double get(int i, int j) override;

	double diag(int i) override;

	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	void Create(double* pv, int* pp, int N);

protected:
	double*	m_pd;			//!< Matrix values
	int*	m_ppointers;	//!< arrays of indices to diagonal elements
};
