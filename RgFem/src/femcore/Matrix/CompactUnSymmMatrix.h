/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "CompactMatrix.h"
#include "femcore/fem_export.h"

struct MatrixItem
{
	int	row, col;
	double	val;
};

//=============================================================================
//! This class stores a general, sparse Matrix in Compact Row Storage format
class FEM_EXPORT CRSSparseMatrix : public CompactMatrix
{
public:
	class Iterator
	{
	public:
		Iterator(CRSSparseMatrix* A);
		bool valid();
		void next();
		void reset();

		MatrixItem get();
		void set(double v);

	private:
		int	r, n;
		CRSSparseMatrix*	m_A;
	};

public:
	//! constructor
	CRSSparseMatrix(int offset = 0);

	//! copy constructor
	CRSSparseMatrix(const CRSSparseMatrix& A);

	//! Create the Matrix structure from the SparseMatrixProfile
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble the element Matrix into the global Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lm) override;

	//! assemble a Matrix into the sparse Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	//! add a value to the Matrix item
	void add(int i, int j, double v) override;

	//! set the Matrix item
	void set(int i, int j, double v) override;

	//! get a Matrix item
	double get(int i, int j) override;

	//! return the diagonal value
	double diag(int i) override;

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

	//! see if a Matrix element is defined
	bool check(int i, int j) override;

	// scale Matrix 
	void scale(double s);
	void scale(const std::vector<double>& L, const std::vector<double>& R) override;

	//! extract a block of this Matrix
	void get(int i0, int j0, int nr, int nc, CSRMatrix& M);

	//! is the Matrix symmetric or not
	bool isSymmetric() override { return false; }

	//! is this a row-based format or not
	bool isRowBased() override { return true; }

	//! calculate the inf norm
	double infNorm() const override;

	//! calculate the one norm
	double oneNorm() const override;

	//! make the Matrix a unit Matrix (retains sparsity pattern)
	void makeUnit();

	//! Create a copy of the Matrix (does not copy values)
	CRSSparseMatrix* Copy(int offset);

	//! Copy the values from another Matrix
	void CopyValues(CompactMatrix* A);

	//! convert to another format (currently only offset can be changed)
	bool Convert(int newOffset);
};

//=============================================================================
//! This class stores a sparse Matrix in Compact Column Storage format

class FEM_EXPORT CCSSparseMatrix : public CompactMatrix
{
public:
	//! constructor
	CCSSparseMatrix(int offset = 0);

	//! copy constructor
	CCSSparseMatrix(const CCSSparseMatrix& A);

	//! Create the Matrix structure from the SparseMatrixProfile
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble the element Matrix into the global Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lm) override;

	//! assemble a Matrix into the sparse Matrix
	void Assemble(const Matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	//! add a value to the Matrix item
	void add(int i, int j, double v) override;

	//! set the Matrix item
	void set(int i, int j, double v) override;

	//! get a Matrix item
	double get(int i, int j) override;

	//! return the diagonal value
	double diag(int i) override;

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

	//! see if a Matrix element is defined
	bool check(int i, int j) override;

	//! is the Matrix symmetric or not
	bool isSymmetric() override { return false; }

	//! is this a row-based format or not
	bool isRowBased() override { return false; }

	//! calculate the inf norm
	double infNorm() const override;

	//! calculate the one norm
	double oneNorm() const override;

	//! do row (L) and column (R) scaling
	void scale(const std::vector<double>& L, const std::vector<double>& R) override;
};
