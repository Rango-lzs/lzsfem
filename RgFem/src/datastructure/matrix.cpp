#include "Matrix.h"
#include "datastructure/Vector3d.h"
#include "datastructure//Matrix3d.h"
#include <assert.h>
#include <math.h>
#include <memory>
using namespace std;

//-----------------------------------------------------------------------------
// These functions are defined in colsol.cpp
void lubksb(double**a, int n, int *indx, double b[]);
void ludcmp(double**a, int n, int* indx);

//-----------------------------------------------------------------------------
// These functions are defined in svd.cpp
void svbksb(Matrix& u, vector<double>& w, Matrix& v, vector<double>& b, vector<double>& x);
void svdcmp(Matrix& a, vector<double>& w, Matrix& v);

//-----------------------------------------------------------------------------
vector<double> operator / (vector<double>& b, Matrix& m)
{
	int n = (int)b.size();

	vector<double> x(b);
	vector<int> indx(n);
	Matrix a(m);

	ludcmp(a, n, &indx[0]);
	lubksb(a, n, &indx[0], &x[0]);

	return x;
}

vector<double> operator * (Matrix& m, vector<double>& b)
{
	int i, j;
	int NR = m.rows();
	int NC = m.columns();
	assert(NC == b.size());
	vector<double> r(NR);
	for (i=0; i<NR; ++i)
	{
		r[i] = 0.0;
		for (j=0; j<NC; ++j) r[i] += m[i][j]*b[j];
	}

	return r;
}

//-----------------------------------------------------------------------------
void Matrix::alloc(int nr, int nc)
{
	m_nr = nr;
	m_nc = nc;
	m_nsize = nr*nc;

	m_pd = new double [m_nsize];
	m_pr = new double*[nr];
	for (int i=0; i<nr; i++) m_pr[i] = m_pd + i*nc;
}

//-----------------------------------------------------------------------------
//! Constructor for Matrix class. 
Matrix::Matrix(int nr, int nc)
{
	alloc(nr, nc);
}

//-----------------------------------------------------------------------------
//! Matrix destructor
void Matrix::clear()
{
	if (m_pr) delete [] m_pr;
	if (m_pd) delete [] m_pd;
	m_pd = 0;
	m_pr = 0;
	m_nr = m_nc = 0;
}

//-----------------------------------------------------------------------------
//! Copy constructor for Matrix class. 
Matrix::Matrix(const Matrix& m)
{
	alloc(m.m_nr, m.m_nc);
	for (int i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];
}

//-----------------------------------------------------------------------------
//! constructor
Matrix::Matrix(const Matrix3d& m)
{
	alloc(3, 3);
	m_pr[0][0] = m[0][0]; m_pr[0][1] = m[0][1]; m_pr[0][2] = m[0][2];
	m_pr[1][0] = m[1][0]; m_pr[1][1] = m[1][1]; m_pr[1][2] = m[1][2];
	m_pr[2][0] = m[2][0]; m_pr[2][1] = m[2][1]; m_pr[2][2] = m[2][2];
}

//-----------------------------------------------------------------------------
//! assignment operator
Matrix& Matrix::operator = (const Matrix3d& m)
{
	resize(3, 3);
	m_pr[0][0] = m[0][0]; m_pr[0][1] = m[0][1]; m_pr[0][2] = m[0][2];
	m_pr[1][0] = m[1][0]; m_pr[1][1] = m[1][1]; m_pr[1][2] = m[1][2];
	m_pr[2][0] = m[2][0]; m_pr[2][1] = m[2][1]; m_pr[2][2] = m[2][2];
	return *this;
}

//-----------------------------------------------------------------------------
Matrix& Matrix::operator = (const Matrix& m)
{
	if ((m.m_nr != m_nr) || (m.m_nc != m_nc))
	{
		clear();
		alloc(m.m_nr, m.m_nc);
	}
	for (int i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];

	return (*this);
}

//-----------------------------------------------------------------------------
void Matrix::resize(int nr, int nc)
{
	if ((nr != m_nr) || (nc != m_nc))
	{
		clear();
		alloc(nr, nc);
	}
}

//-----------------------------------------------------------------------------
Matrix Matrix::operator * (double a) const
{
	Matrix m(m_nr, m_nc);
	int n = m_nr*m_nc;
	for (int i = 0; i < n; ++i) m.m_pd[i] = m_pd[i] * a;
	return m;
}

//-----------------------------------------------------------------------------
Matrix Matrix::operator * (const Matrix& m) const
{
	assert(m_nc == m.m_nr);
	Matrix a(m_nr, m.m_nc);

	// NOTE: commented this out since this can be called for small matrices.
	// TODO: make a separate function that implements a parallel version. (e.g. pmult) 
//	#pragma omp parallel for shared(a)
	for (int i = 0; i < m_nr; ++i)
	{
		double* pa = a.m_pr[i];
		for (int j = 0; j < m.m_nc; ++j) pa[j] = 0.0;
		for (int k = 0; k < m_nc; ++k)
		{
			const double pik = m_pr[i][k];
			const double* pm = m.m_pr[k];
			for (int j = 0; j < m.m_nc; ++j)
			{
				pa[j] += pik * pm[j];
			}
		}
	}
	return a;
}

//-----------------------------------------------------------------------------
// Calculate the LU decomposition of this Matrix. Note that this will modify
// the Matrix. This is used for repeated solves of a linear system. Use 
// lusolve for solving after lufactor.
void Matrix::lufactor(vector<int>& indx)
{
	// make sure this is a square Matrix
	assert(m_nr == m_nc);

	// do a LU decomposition
	int n = m_nr;
	indx.resize(n);
	ludcmp(*(this), n, &indx[0]);
}

//-----------------------------------------------------------------------------
// Solve the linear system Ax=b, where A has been factored using lufactor.
// The indx array is the same one that was returned from lufactor
void Matrix::lusolve(vector<double>& b, vector<int>& indx)
{
	// make sure this is a square Matrix
	assert(m_nr == m_nc);
	lubksb(*(this), m_nr, &indx[0], &b[0]);
}

//-----------------------------------------------------------------------------
Matrix Matrix::inverse()
{
	// make sure this is a square Matrix
	assert(m_nr == m_nc);

	// make a copy of this Matrix
	// since we don't want to change it
	Matrix a(*this);

	// do a LU decomposition
	int n = m_nr;
	vector<int> indx(n);
	ludcmp(a, n, &indx[0]);

	// allocate the inverse Matrix
	Matrix ai(n, n);

	// do a backsubstituation on the columns of a
	vector<double> b; b.assign(n, 0);
	for (int j=0; j<n; ++j)
	{
		b[j] = 1;
		lubksb(a, n, &indx[0], &b[0]);

		for (int i=0; i<n; ++i)
		{
			ai[i][j] = b[i];
			b[i] = 0;
		}
	}

	return ai;
}

//-----------------------------------------------------------------------------
// Matrix using singular value decomposition
Matrix Matrix::svd_inverse()
{
	Matrix U(*this);
	Matrix V(m_nr, m_nc);
	vector<double> w(m_nc);

	// calculate the decomposition
	svdcmp(U, w, V);

	Matrix Ai(m_nc, m_nr); // inverse
	for (int i=0; i<m_nc; ++i)
		for (int j=0; j<m_nr; ++j)
		{
			double s = 0.0;
			for (int k=0;k<m_nc; ++k)
			{
				if (w[k] > 0.0)
				{
					s += (V[i][k]*U[j][k])/w[k];
				}
			}
			Ai[i][j] = s;
		}
	return Ai;
}

//-----------------------------------------------------------------------------
Matrix Matrix::transpose()
{
	int i, j;
	Matrix At(m_nc, m_nr);
	for (i=0; i<m_nr; ++i)
		for (j=0; j<m_nc; ++j) At[j][i] = m_pr[i][j];
	return At;
}

//-----------------------------------------------------------------------------
Matrix& Matrix::operator += (const Matrix& m)
{
	assert((m_nr == m.m_nr ) && (m_nc == m.m_nc));
	for (int i=0; i<m_nsize; ++i) m_pd[i] += m.m_pd[i];
	return (*this);
}

//-----------------------------------------------------------------------------
Matrix& Matrix::operator -= (const Matrix& m)
{
	assert((m_nr == m.m_nr ) && (m_nc == m.m_nc));
	for (int i=0; i<m_nsize; ++i) m_pd[i] -= m.m_pd[i];
	return (*this);
}

//-----------------------------------------------------------------------------
Matrix Matrix::operator * (const Vector3d& v) const
{
	assert(m_nc == 3);
	Matrix A(m_nr, 1);
	const Matrix& T = *this;
	for (int i = 0; i < m_nr; ++i)
	{
		A[i][0] = T[i][0]*v.x + T[i][1] * v.y + T[i][2] * v.z;
	}
	return A;
}

//-----------------------------------------------------------------------------
// calculate outer product of a vector to produce a Matrix
Matrix outer_product(vector<double>& a)
{
	int n = (int) a.size();
	Matrix m(n,n);
	for (int i=0; i<n; ++i)
	{
		for (int j=0; j<n; ++j) m[i][j] = a[i]*a[j];
	}
	return m;
}

//-----------------------------------------------------------------------------
// Calculates the infinity norm. That is, the max of the absolute row sum.
double Matrix::inf_norm()
{
	Matrix& self = (*this);
	double m = 0;
	for (int j=0; j<m_nr; ++j)
	{
		double s = 0;
		for (int i=0; i<m_nc; ++i) s += fabs(self[j][i]);
		if (s > m) m = s;
	}
	return m;
}

//-----------------------------------------------------------------------------
void Matrix::get(int i, int j, int rows, int cols, Matrix& A) const
{
	// make sure we create a valid Matrix
	if ((rows <= 0) || (cols <= 0)) return;

	// initialize the Matrix
	A.resize(rows, cols);
	A.zero();

	// make sure the bounds are within this Matrix
	if ((i >= m_nr) || (j >= m_nc)) return;
	if ((i + rows <= 0) || (j + cols <= 0)) return;

	// set range
	int r0 = 0;
	int r1 = rows - 1;
	int c0 = 0;
	int c1 = cols - 1;
	if (i < 0) r0 -= i;
	if (j < 0) c0 -= j;
	if (i + rows > m_nr) r1 -= rows + i - m_nr;
	if (j + cols > m_nc) c1 -= cols + j - m_nc;

	for (int r=r0; r<=r1; ++r)
		for (int c=c0; c<=c1; ++c)
				A[r][c] = (*this)(i+r, j+c);
}

//-----------------------------------------------------------------------------
// fill a Matrix
void Matrix::fill(int i, int j, int rows, int cols, double val)
{
	if ((i >= m_nr) || (j >= m_nc)) return;
	if ((i + rows <= 0) || (j + cols <= 0)) return;

	int r0 = 0;
	int r1 = rows - 1;
	int c0 = 0;
	int c1 = cols - 1;
	if (i < 0) r0 -= i;
	if (j < 0) c0 -= j;
	if (i + rows > m_nr) r1 -= rows + i - m_nr;
	if (j + cols > m_nc) c1 -= cols + j - m_nc;

	for (int r = r0; r<=r1; ++r)
		for (int c = c0; c<=c1; ++c)
			(*this)(i + r, j + c) = val;
}

//-----------------------------------------------------------------------------
// solve the linear system Ax=b
void Matrix::solve(vector<double>& x, const vector<double>& b)
{
	Matrix A(*this);
	vector<int> index;
	A.lufactor(index);
	x = b;
	A.lusolve(x, index);
}

//-----------------------------------------------------------------------------
Matrix Matrix::operator + (const Matrix& m) const
{
	Matrix s(*this);
	for (int i=0; i<m_nr; ++i)
		for (int j=0; j<m_nc; ++j)
			s(i,j) += m(i,j);

	return s;
}

//-----------------------------------------------------------------------------
Matrix Matrix::operator - (const Matrix& m) const
{
	Matrix s(*this);
	for (int i = 0; i<m_nr; ++i)
		for (int j = 0; j<m_nc; ++j)
			s(i, j) -= m(i, j);

	return s;
}

//-----------------------------------------------------------------------------
void Matrix::mult(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < m_nr; ++i)
	{
		double* di = m_pd + i * m_nc;
		y[i] = 0.0;
		for (int j = 0; j < m_nc; ++j) y[i] += di[j] * x[j];
	}
}

// Here y = this*m*x. This is useful if this*m is a very large Matrix, 
// but is then immediately multiplied by a vector, brining its size down 
// to m_nr x 1. In this unique circumstance, the memory requrements can be 
// drastically lower. 
void Matrix::mult(const Matrix& m, std::vector<double>& x, std::vector<double>& y)
{
    assert(m_nc == m.m_nr);
    assert(m_nr == x.size());
    
    #pragma omp parallel for
	for (int i = 0; i < m_nr; ++i)
	{
        vector<double> temp(m.m_nc, 0.0);
		for (int k = 0; k < m_nc; ++k)
		{
			const double pik = m_pr[i][k];
			const double* pm = m.m_pr[k];
			for (int j = 0; j < m.m_nc; ++j)
			{
				temp[j] += pik * pm[j];
			}
		}

        y[i] = 0.0;
        for(int j = 0; j < m_nr; j++)
        {
            y[i] += temp[j] * x[j];
        }
    
	}

}

//-----------------------------------------------------------------------------
void Matrix::mult_transpose(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < m_nc; ++i) y[i] = 0.0;

	for (int i = 0; i < m_nr; ++i)
	{
		double* di = m_pd + i * m_nc;
		for (int j = 0; j < m_nc; ++j) y[j] += di[j] * x[i];
	}
}

//-----------------------------------------------------------------------------
void Matrix::mult_transpose_self(Matrix& AAt)
{
	Matrix& A = *this;
	int N = m_nc;
	int R = m_nr;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			double& aij = AAt[i][j];
			aij = 0.0;
			for (int k = 0; k < R; ++k) aij += A[k][i] * A[k][j];
		}
}

//-----------------------------------------------------------------------------
bool Matrix::lsq_solve(vector<double>& x, vector<double>& b)
{
	if ((int)x.size() != m_nc) return false;
	if ((int)b.size() != m_nr) return false;

	vector<double> y(m_nc);
	mult_transpose(b, y);

	Matrix AA(m_nc, m_nc);
	mult_transpose_self(AA);

	AA.solve(x, y);

	return true;
}

#define ROTATE(a, i, j, k, l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l] = h + s*(g - h*tau);

//-----------------------------------------------------------------------------
bool Matrix::eigen_vectors(Matrix& Eigen, vector<double>& eigen_values)
{
	Matrix& A = *this;
	int N = m_nc;
	int R = m_nr;
	const int NMAX = 50;
	double sm, tresh, g, h, t, c, tau, s, th;
	const double eps = 0;//1.0e-15;
	int k;

	//initialize Eigen to identity
	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < N; j++) Eigen[i][j] = 0;
		Eigen[i][i] = 1;
	}
	vector<double> b;
	b.reserve(R);
	vector<double> z;
	z.reserve(R);

	// initialize b and eigen_values to the diagonal of A
	for (int i = 0; i < R; i++)
	{
		b.push_back(A[i][i]);
		eigen_values.push_back(A[i][i]);
		z.push_back(0);
	}
	// loop
	int n, nrot = 0;
	for (n = 0; n < NMAX; ++n)
	{
		// sum off-diagonal elements
		sm = 0;
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++)
				sm += fabs(A[i][j]);
		}
		if (sm <= eps) {
			break;
		}
		// set the treshold
		if (n < 3) tresh = 0.2 * sm / (R * R); else tresh = 0.0;

		// loop over off-diagonal elements
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {

				g = 100.0 * fabs(A[i][j]);
				// after four sweeps, skip the rotation if the off-diagonal element is small
				if ((n > 3) && ((fabs(eigen_values[i]) + g) == fabs(eigen_values[i])) && ((fabs(eigen_values[j]) + g) == fabs(eigen_values[j])))
				{
					A[i][j] = 0.0;
				}
				else if (fabs(A[i][j]) > tresh) {
					h = eigen_values[j] - eigen_values[i];
					if ((fabs(h) + g) == fabs(h))
						t = A[i][j] / h;
					else
					{
						th = 0.5 * h / A[i][j];
						t = 1.0 / (fabs(th) + sqrt(1 + th * th));
						if (th < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1.0 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * A[i][j];
					z[i] -= h;
					z[j] += h;
					eigen_values[i] -= h;
					eigen_values[j] += h;
					A[i][j] = 0;
					for (k = 0; k <= i - 1; ++k) { ROTATE(A, k, i, k, j) }
					for (k = i + 1; k <= j - 1; ++k) { ROTATE(A, i, k, k, j) }
					for (k = j + 1; k < N; ++k) { ROTATE(A, i, k, j, k) }
					for (k = 0; k < N; ++k) { ROTATE(Eigen, k, i, k, j) }
					++nrot;
				}
			}
		}//end of for loop

		//Update eigen_values with the sum. Reinitialize z.
		for (int i = 0; i < R; ++i)
		{
			b[i] += z[i];
			eigen_values[i] = b[i];
			z[i] = 0.0;
		}
	}

	// we sure we converged
	assert(n < NMAX);
	return true;
}
