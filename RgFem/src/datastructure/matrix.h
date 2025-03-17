#pragma once

#include "femcore/fem_export.h"
#include <memory.h>
#include <vector>

class Matrix3d;
class Vector3d;

//-----------------------------------------------------------------------------
//!Ò»°ãµÄ¾ØÕó£¬
class FEM_EXPORT Matrix
{
public:
	//! constructor
	Matrix() : m_nr(0), m_nc(0), m_nsize(0), m_pd(nullptr), m_pr(nullptr) {}

	//! constructor
	Matrix(int nr, int nc);

	//! copy constructor
	Matrix(const Matrix& m);

	//! constructor
	Matrix(const Matrix3d& m);

	//! move constructor
	Matrix(Matrix&& m);

	//! assignment operator
	Matrix& operator = (const Matrix& m);

	//! move assigment operator
	Matrix& operator = (Matrix&& m);

	//! assignment operator
	Matrix& operator = (const Matrix3d& m);

	//! Matrix reallocation
	void resize(int nr, int nc);

	//! destructor
	~Matrix() { clear(); }

	//! access operator
	double * operator [] (int l) { return m_pr[l]; }
	const double* operator [] (int l) const { return m_pr[l]; }
	double& operator () (int i, int j) { return m_pr[i][j]; }
	double operator () (int i, int j) const { return m_pr[i][j]; }
	operator double** () { return m_pr; }

	int rows   () const { return m_nr; }
	int columns() const { return m_nc; }

	void zero() { memset(m_pd, 0, sizeof(double)*m_nsize); }

	//! Matrix transpose
	Matrix transpose();

	//! Matrix inversion
	Matrix inverse();

	//! Matrix inverse using SVD
	Matrix svd_inverse();

	//! Matrix operators
	Matrix operator *(double a) const;
	Matrix operator * (const Matrix& m) const;

	Matrix operator + (const Matrix& m) const;

	Matrix operator - (const Matrix& m) const;

	Matrix& operator += (const Matrix& m);

	Matrix& operator -= (const Matrix& m);

	Matrix& operator *=(double g)
	{
		for (int i=0; i<m_nsize; ++i) m_pd[i] *= g;
		return *this;
	}

	Matrix operator * (const Vector3d& v) const;

	// calculate the LU decomposition
	// note that this modifies the Matrix
	void lufactor(std::vector<int>& indx);

	// solve using the lu factor calculated with lufactor
	void lusolve(std::vector<double>& b, std::vector<int>& indx);

	// solve the linear system Ax=b
	void solve(std::vector<double>& x, const std::vector<double>& b);

	bool lsq_solve(std::vector<double>& x, std::vector<double>& b);

	bool eigen_vectors(Matrix& Eigen, std::vector<double>& eigen_values);

	// infinity-norm
	double inf_norm();

	void mult(std::vector<double>& x, std::vector<double>& y);
    void mult(const Matrix& m, std::vector<double>& x, std::vector<double>& y);
	void mult_transpose(std::vector<double>& x, std::vector<double>& y);
	void mult_transpose_self(Matrix& AAt);


	// extract a Matrix block
	// the returned Matrix will have the dimensions rows x cols
	// if the Matrix doesn't fit in this Matrix, the missing entries will be set to zero
	void get(int i, int j, int rows, int cols, Matrix& A) const;

	// fill a Matrix
	void fill(int i, int j, int rows, int cols, double val);

private:
	void alloc(int nr, int nc);
	void clear();

protected:
	double**	m_pr;	// pointer to rows
	double*		m_pd;	// Matrix elements

	int	m_nr;		// nr of rows
	int	m_nc;		// nr of columns
	int	m_nsize;	// size of Matrix (ie. total nr of elements = nr*nc)
};

std::vector<double> FEM_EXPORT operator / (std::vector<double>& b, Matrix& m);
std::vector<double> FEM_EXPORT operator * (Matrix& m, std::vector<double>& b);
Matrix FEM_EXPORT outer_product(std::vector<double>& a);

//! move constructor
inline Matrix::Matrix(Matrix&& m)
{
	m_nr = m.m_nr;
	m_nc = m.m_nc;
	m_pd = m.m_pd;
	m_pr = m.m_pr;

	m.m_pr = nullptr;
	m.m_pd = nullptr;
}

//! move assigment operator
inline Matrix& Matrix::operator = (Matrix&& m)
{
	if (this != &m)
	{
		if (m_pd) delete[] m_pd;
		if (m_pr) delete[] m_pr;

		m_nr = m.m_nr;
		m_nc = m.m_nc;
		m_pd = m.m_pd;
		m_pr = m.m_pr;

		m.m_pr = nullptr;
		m.m_pd = nullptr;
	}

	return *this;
}
