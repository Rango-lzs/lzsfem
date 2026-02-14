#pragma once

#include "SparseMatrix.h"
#include "femcore/NewtonSolver/SolverBase.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FEMesh;
class FESurface;
class RgElement;

//-----------------------------------------------------------------------------
//! This class represents an element Matrix, i.e. a Matrix of values and the row and
//! column indices of the corresponding Matrix elements in the global Matrix. 
class FEM_EXPORT FEElementMatrix : public Matrix
{
public:
	// default constructor
	FEElementMatrix(){}
	FEElementMatrix(int nr, int nc) : Matrix(nr, nc) {}
	FEElementMatrix(const RgElement& el);

	// constructor for symmetric matrices
	FEElementMatrix(const RgElement& el, const std::vector<int>& lmi);

	// constructor
	FEElementMatrix(const RgElement& el, std::vector<int>& lmi, std::vector<int>& lmj);

	// copy constructor
	FEElementMatrix(const FEElementMatrix& ke);
	FEElementMatrix(const FEElementMatrix& ke, double scale);

	// assignment operator
	void operator = (const Matrix& ke);

	// row indices
	std::vector<int>& RowIndices() { return m_lmi; }
	const std::vector<int>& RowIndices() const { return m_lmi; }

	// column indices
	std::vector<int>& ColumnsIndices() { return m_lmj; }
	const std::vector<int>& ColumnsIndices() const { return m_lmj; }

	// set the row and columnd indices (assuming they are the same)
	void SetIndices(const std::vector<int>& lm) { m_lmi = m_lmj = lm; }

	// set the row and columnd indices
	void SetIndices(const std::vector<int>& lmr, const std::vector<int>& lmc) { m_lmi = lmr; m_lmj = lmc; }

	// Set the node indices
	void SetNodes(const std::vector<int>& en) { m_node = en; }

	// get the nodes
	const std::vector<int>& Nodes() const { return m_node; }

private:
	std::vector<int>	m_node;	//!< node indices
	std::vector<int>	m_lmi;	//!< row indices
	std::vector<int>	m_lmj;	//!< column indices
};

//-----------------------------------------------------------------------------
//! This class implements a global system Matrix.

//! The global system Matrix is usually created by the discretization of the FE 
//! equations into a linear system of equations. The structure of it depends greatly
//! on the element connectivity and usually results in a sparse Matrix structure. 
//! Several sparse Matrix structures are supported (Compact, Skyline, etc.) and to 
//! simplify the creation of the specific Matrix structure, the FEGlobalMatrix offers
//! functionality to create the global Matrix structure without the need to know 
//! what particular sparse Matrix format is used by the linear solver.

//! \todo I think the SparseMatrixProfile can handle all of the build functions.

class FEM_EXPORT FEGlobalMatrix
{
protected:
	enum { MAX_LM_SIZE = 64000 };

public:
	//! constructor
	FEGlobalMatrix(SparseMatrix* pK, bool del = true);

	//! destructor
	virtual ~FEGlobalMatrix();

	//! construct the stiffness Matrix from a FEM object
	bool Create(FEModel* pfem, int neq, bool breset);

	//! construct the stiffness Matrix from a mesh
	bool Create(FEMesh& mesh, int neq);

	//! construct the stiffness Matrix from a mesh
	bool Create(FEMesh& mesh, int nstart, int nend);

	//! construct a stiffness Matrix from a surface
	//! The boundary array is a list of equation numbers.
	bool Create(const FESurface& surf, const std::vector<int>& equationIDs);

	//! clears the sparse Matrix that stores the stiffness Matrix
	void Clear();

	//! Assembly routine
	virtual void Assemble(const FEElementMatrix& ke);

	//! return the nonzeroes in the sparse Matrix
	int NonZeroes() { return m_pA->NonZeroes(); }

	//! return the number of rows
	int Rows() { return m_pA->Rows(); }

	//! converts a FEGlobalMatrix to a SparseMatrix
	operator SparseMatrix* () { return m_pA; }

	//! converts a FEGlobalMatrix to a SparseMatrix
	operator SparseMatrix& () { return *m_pA;}

	//! return a pointer to the sparse Matrix
	SparseMatrix* GetSparseMatrixPtr() { return m_pA; }

	//! zero the sparse Matrix
	void Zero() { m_pA->Zero(); }

	//! get the sparse Matrix profile
	SparseMatrixProfile* GetSparseMatrixProfile() { return m_pMP; }

public:
	void build_begin(int neq);
	void build_add(std::vector<int>& lm);
	void build_end();
	void build_flush();

protected:
	SparseMatrix*	m_pA;	//!< the actual global stiffness Matrix
	bool			m_delA;	//!< delete A in destructor

	// The following data structures are used to incrementally
	// build the profile of the sparse Matrix

	SparseMatrixProfile*	m_pMP;		//!< profile of sparse Matrix
	SparseMatrixProfile		m_MPs;		//!< the "static" part of the Matrix profile
	std::vector< std::vector<int> >	m_LM;		//!< used for building the stiffness Matrix
	int	m_nlm;				//!< nr of elements in m_LM array
};
