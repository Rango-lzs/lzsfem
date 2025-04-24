#pragma once
#include "femcore/Matrix/MatrixProfile.h"
#include "datastructure/MatrixOperator.h"
#include "datastructure/Matrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! Base class for sparse matrices

//! This is the base class for the sparse Matrix classes and defines the interface
//! to the different Matrix classes

class FEM_EXPORT SparseMatrix : public MatrixOperator
{
public:
	//! constructor
	SparseMatrix();

	//! destructor
	virtual ~SparseMatrix();

public:
	//! return number of rows
	int Rows() const { return m_nrow; }

	//! return number of columns
	int Columns() const { return m_ncol; }

	//! is the Matrix square?
	bool IsSquare() const { return (m_nrow == m_ncol); }

	//! return number of nonzeros
	int NonZeroes() const { return m_nsize; }

public: // functions to be overwritten in derived classes

	//! set all Matrix elements to zero
	virtual void Zero() = 0;

	//! Create a sparse Matrix from a sparse-Matrix profile
	virtual void Create(SparseMatrixProfile& MP) = 0;

	//! assemble a Matrix into the sparse Matrix
	virtual void Assemble(const Matrix& ke, const std::vector<int>& lm) = 0;

	//! assemble a Matrix into the sparse Matrix
	virtual void Assemble(const Matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) = 0;

	//! check if an entry was allocated
	virtual bool check(int i, int j) = 0;

	//! set entry to value
	virtual void set(int i, int j, double v) = 0;

	//! add value to entry
	virtual void add(int i, int j, double v) = 0;

	//! retrieve value
	virtual double get(int i, int j) { return 0; }

	//! get the diagonal value
	virtual double diag(int i) = 0;

	//! release memory for storing data
	virtual void Clear();

	//! scale Matrix
	virtual void scale(const std::vector<double>& L, const std::vector<double>& R);

public:
	//! multiply with vector
	bool mult_vector(double* x, double* r) override { assert(false); return false; }

public:
	// NOTE: The following functions are only used by the compact matrices, but I need to be able to override them
	// for the JFNKMatrix so I've moved them here. 
	virtual double* Values() { return 0; }
	virtual int*    Indices() { return 0; }
	virtual int*    Pointers() { return 0; }
	virtual int     Offset() const { return 0; }

protected:
	// NOTE: These values are set by derived classes
	int	m_nrow, m_ncol;		//!< dimension of Matrix
	int	m_nsize;			//!< number of nonzeroes (i.e. Matrix elements actually allocated)
};
