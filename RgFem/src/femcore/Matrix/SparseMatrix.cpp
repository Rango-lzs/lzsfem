#include "SparseMatrix.h"
#include <memory.h>
using namespace std;

//-----------------------------------------------------------------------------
SparseMatrix::SparseMatrix()
{
	m_nrow = m_ncol = 0;
	m_nsize = 0;
}

SparseMatrix::~SparseMatrix()
{
}

void SparseMatrix::Clear()
{
	m_nrow = m_ncol = 0;
	m_nsize = 0;
}

//! scale matrix
void SparseMatrix::scale(const vector<double>& L, const vector<double>& R)
{
	assert(false);
}
