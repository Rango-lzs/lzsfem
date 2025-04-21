#pragma once
#include <vector>

// template class for storing 2D data
template <typename T> class Table2d
{
public:
	// constructors
	Table2d() : m_rows(0), m_cols(0) {}
	Table2d(int nrows, int ncols) : m_rows(0), m_cols(0) { resize(nrows, ncols); }
	Table2d(const Table2d& t) { m_data = t.m_data; m_rows = t.m_rows; m_cols = t.m_cols; }

	// assignment operator
	Table2d& operator = (const Table2d& t) { m_data = t.m_data; m_rows = t.m_rows; m_cols = t.m_cols; return (*this); }

	// resize Table2d
	void resize(int nrows, int ncols, T def = T(0))
	{
		std::vector<T> tmp(m_data);
		m_data.assign(nrows*ncols, def);

		int nr = (nrows < m_rows ? nrows : m_rows);
		int nc = (ncols < m_cols ? ncols : m_cols);
		for (int i=0; i<nr; ++i)
		{
			for (int j=0; j<nc; ++j) m_data[i*ncols + j] = tmp[i*m_cols + j];
		}

		m_rows = nrows;
		m_cols = ncols;
	}

	// assign a value to the entire Table2d
	void set(const T& v) 
	{ 
		if (m_data.empty() == false) 
			m_data.assign(m_rows*m_cols, v); 
	}

	// get sizes
	int rows() const { return m_rows; }
	int columns() const { return m_cols; }

	// access operator
	const T& operator () (int i, int j) const { return m_data[i*m_cols + j]; }
	T& operator () (int i, int j) { return m_data[i*m_cols + j]; }

private:
	std::vector<T>	m_data;
	int				m_rows, m_cols;
};
