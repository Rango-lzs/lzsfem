#include "CSRMatrix.h"
#include <assert.h>

CSRMatrix::CSRMatrix() : m_nr(0), m_nc(0), m_offset(0)
{
}

// create a Matrix of given size
CSRMatrix::CSRMatrix(int rows, int cols, int noffset) : m_nr(rows), m_nc(cols), m_offset(noffset)
{
	m_rowIndex.resize(rows+1, m_offset);
}

// Create Matrix
void CSRMatrix::create(int nr, int nc, int noffset)
{
	m_nr = nr;
	m_nc = nc;
	m_offset = noffset;
	m_rowIndex.resize(nr+1, noffset);
	m_columns.clear();
	m_values.clear();
}

// copy constructor
CSRMatrix::CSRMatrix(const CSRMatrix& A)
{
	m_nr = A.m_nr;
	m_nc = A.m_nc;
	m_offset = A.m_offset;
	m_rowIndex = A.m_rowIndex;
	m_columns = A.m_columns;
	m_values = A.m_values;
}

// assignment operator
void CSRMatrix::operator = (const CSRMatrix& A)
{
	m_nr = A.m_nr;
	m_nc = A.m_nc;
	m_offset = A.m_offset;
	m_rowIndex = A.m_rowIndex;
	m_columns = A.m_columns;
	m_values = A.m_values;
}

// sets a value, inserting it if necessary
void CSRMatrix::set(int i, int j, double val)
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	// get the start column index for the given row and the non-zero count for that row
	int col = m_rowIndex[i] - m_offset;
	int count = m_rowIndex[i+1] - m_rowIndex[i];

	// see if the row is empty
	if (count == 0)
	{
		m_columns.insert(m_columns.begin() + col, j);
		m_values.insert(m_values.begin() + col, val);
	}
	else
	{
		// see if this column would be before the first entry
		if (j < m_columns[col] - m_offset)
		{
			m_columns.insert(m_columns.begin() + col, j);
			m_values.insert(m_values.begin() + col, val);
		}
		// see if this column would be past the last entry
		else if (j > m_columns[col + count - 1] - m_offset)
		{
			m_columns.insert(m_columns.begin() + col + count, j);
			m_values.insert(m_values.begin() + col + count, val);
		}
		else {
			// find the column index
			for (int n=0; n<count; ++n)
			{
				// see if it alreay exists
				if (m_columns[col + n] - m_offset == j)
				{
					// if so, replace the value and return
					m_values[col+n] = val;
					return;
				}
				else if (m_columns[col + n] - m_offset > j)
				{
					// we found an index that is bigger so insert this value before this one
					m_columns.insert(m_columns.begin() + col + n, j);
					m_values.insert(m_values.begin() + col + n, val);
					break;
				}
			}
		}
	}

	// increase row counts
	for (int n = i + 1; n <= m_nr; ++n) m_rowIndex[n]++;
	assert(m_rowIndex[m_nr] == m_values.size());
	assert(m_columns.size() == m_values.size());
}

// get a value
double CSRMatrix::operator () (int i, int j) const
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	int col = m_rowIndex[i] - m_offset;
	int count = m_rowIndex[i + 1] - m_rowIndex[i];

	if (count == 0) return 0.0;
	if (j < m_columns[col] - m_offset) return 0.0;
	if (j > m_columns[col + count - 1] - m_offset) return 0.0;
	for (int n=0; n<count; ++n)
	{
		if (m_columns[col + n] - m_offset == j) return m_values[col + n];
	}
	return 0.0;
}

// see if a Matrix entry was allocated
bool CSRMatrix::isAlloc(int i, int j) const
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	int col = m_rowIndex[i] - m_offset;
	int count = m_rowIndex[i + 1] - m_rowIndex[i];

	if (count == 0) return false;
	if (j < m_columns[col] - m_offset) return false;
	if (j > m_columns[col + count - 1] - m_offset) return false;
	for (int n = 0; n<count; ++n)
	{
		if (m_columns[col + n] - m_offset == j) return true;
	}
	return false;
}

void CSRMatrix::multv(const std::vector<double>& x, std::vector<double>& r)
{
	multv(&x[0], &r[0]);
}

void CSRMatrix::multv(const double* x, double* r)
{
	const int nr = rows();
	const int nc = cols();

	// loop over all rows
	for (int i = 0; i<nr; ++i)
	{
		int col = m_rowIndex[i] - m_offset;
		int count = m_rowIndex[i + 1] - m_rowIndex[i];

		double* pv = &m_values[col];
		int* pi = &m_columns[col];
		double ri = 0.0;
		for (int j = 0; j<count; ++j) ri += pv[j] * x[pi[j] - m_offset];
		r[i] = ri;
	}
}
