#include "LUSolver.h"
#include <math.h>

//-----------------------------------------------------------------------------
LUSolver::LUSolver(FEModel* fem) : LinearSolver(fem), m_pA(nullptr)
{

}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* LUSolver::CreateSparseMatrix(Matrix_Type ntype)
{ 
	return (m_pA = new FECore::DenseMatrix()); 
}

//-----------------------------------------------------------------------------
void LUSolver::SetMatrix(FECore::DenseMatrix* pA) 
{ 
	m_pA = pA; 
}

//-----------------------------------------------------------------------------
bool LUSolver::PreProcess()
{
	// We don't need to do any preprocessing for this solver
	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool LUSolver::Factor()
{
	FECore::DenseMatrix& a = *m_pA;

	const double TINY = 1.0e-20;
	int i, imax, j, k;
	double big, dum, sum, temp;

	int n = a.Rows();
	// create index vector
	indx.resize(n);

	vector<double> vv(n);
	for (i=0; i<n; ++i)
	{
		big = 0;
		for (j=0; j<n; ++j)
			if ((temp=fabs(a(i,j))) > big) big = temp;
		if (big == 0) return false; // singular matrix
		vv[i] = 1.0 / big;
	}

	for (j=0; j<n; ++j)
	{
		for (i=0; i<j; ++i)
		{
			sum = a(i,j);
			for (k=0; k<i; ++k) sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
		}
		big = 0;
		imax = j;
		for (i=j;i<n;++i)
		{
			sum = a(i,j);
			for (k=0; k<j; ++k) sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
			if ((dum=vv[i]*fabs(sum))>=big)
			{
				big = dum;
				imax = i;
			}
		}

		if (j != imax)
		{
			for (k=0; k<n; ++k)
			{
				dum = a(imax,k);
				a(imax,k) = a(j,k);
				a(j,k) = dum;
			}
			vv[imax] = vv[j];
		}

		indx[j] = imax;
		if (a(j,j) == 0) a(j,j) = TINY;
		if (j != n-1)
		{
			dum = 1.0/a(j,j);
			for (i=j+1;i<n; ++i) a(i,j) *= dum;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool LUSolver::BackSolve(double* x, double* b)
{
	FECore::DenseMatrix& a = *m_pA;

	int n = a.Rows();
	for (int i=0; i<n; i++) x[i] = b[i];

	int ii=0;
	for (int i=0; i<n; ++i)
	{
		int ip = indx[i];
		double sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
			for (int j=ii-1;j<i;++j) sum -= a(i,j)*x[j];
		else if (sum != 0)
			ii = i+1;
		x[i] = sum;
	}

	for (int i=n-1; i>=0; --i)
	{
		double sum = x[i];
		for (int j=i+1; j<n; ++j) sum -= a(i,j)*x[j];
		x[i] = sum/a(i,i);
	}

	return false;
}

//-----------------------------------------------------------------------------
void LUSolver::Destroy()
{
	// nothing to destroy
	LinearSolver::Destroy();
}
