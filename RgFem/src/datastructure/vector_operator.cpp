
#include <assert.h>
#include "vector_operator.h"
#include "femcore/FEMesh.h"
#include "femcore/FEDofList.h"
#include <algorithm>

double operator*(const std::vector<double>& a, const std::vector<double>& b)
{
	double sum_p = 0, sum_n = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		double ab = a[i] * b[i];
		if (ab >= 0.0) sum_p += ab; else sum_n += ab;
	}
	return sum_p + sum_n;
}

std::vector<double> operator - (std::vector<double>& a, std::vector<double>& b)
{
	std::vector<double> c(a);
	int n = (int) c.size();
	for (int i=0; i<n; ++i) c[i] -= b[i];
	return c;
}

void operator += (std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}

void operator -= (std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
}

void operator *= (std::vector<double>& a, double b)
{
	for (size_t i=0; i<a.size(); ++i) a[i] *= b;
}

void vcopys(std::vector<double>& a, const std::vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i=0; i<a.size(); ++i) a[i] = b[i]*s;
}

void vadds(std::vector<double>& a, const std::vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] += b[i] * s;
}

void vsubs(std::vector<double>& a, const std::vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] -= b[i] * s;
}

void vscale(std::vector<double>& a, const std::vector<double>& s)
{
	assert(a.size() == s.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] *= s[i];
}

void vsub(std::vector<double>& a, const std::vector<double>& l, const std::vector<double>& r)
{
	assert((a.size()==l.size())&&(a.size()==r.size()));
	for (size_t i=0; i<a.size(); ++i) a[i] = l[i] - r[i];
}

std::vector<double> operator + (const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());
	std::vector<double> s(a);
	for (size_t i = 0; i < s.size(); ++i) s[i] += b[i];
	return s;
}

std::vector<double> operator*(const std::vector<double>& a, double g)
{
	std::vector<double> s(a.size());
	for (size_t i = 0; i < s.size(); ++i) s[i] = a[i]*g;
	return s;
}

std::vector<double> FEM_EXPORT operator - (const std::vector<double>& a)
{
	std::vector<double> s(a.size());
	for (size_t i = 0; i < s.size(); ++i) s[i] = -a[i];
	return s;
}

void gather(std::vector<double>& v, FEMesh& mesh, int ndof)
{
	const int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_dofs[ndof]; if (n >= 0) v[n] = node.get(ndof);
	}
}

void gather(std::vector<double>& v, FEMesh& mesh, const std::vector<int>& dof)
{
	const int NN = mesh.Nodes();
	const int NDOF = (const int) dof.size();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<NDOF; ++j)
		{
			int n = node.m_dofs[dof[j]]; 
			if (n >= 0) v[n] = node.get(dof[j]);
		}
	}
}

void scatter(std::vector<double>& v, FEMesh& mesh, int ndof)
{
	const int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_dofs[ndof];
		if (n >= 0) node.set(ndof, v[n]);
	}
}

void scatter3(std::vector<double>& v, FEMesh& mesh, int ndof1, int ndof2, int ndof3)
{
	const int NN = mesh.Nodes();
#pragma omp parallel for 
	for (int i = 0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n;
		n = node.m_dofs[ndof1]; if (n >= 0) node.set(ndof1, v[n]);
		n = node.m_dofs[ndof2]; if (n >= 0) node.set(ndof2, v[n]);
		n = node.m_dofs[ndof3]; if (n >= 0) node.set(ndof3, v[n]);
	}
}

void scatter(std::vector<double>& v, FEMesh& mesh, const FEDofList& dofs)
{
	const int NN = mesh.Nodes();
	for (int i = 0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < dofs.Size(); ++j)
		{
			int n = node.m_dofs[dofs[j]]; if (n >= 0) node.set(dofs[j], v[n]);
		}
	}
}

double l2_norm(const std::vector<double>& v)
{
	double s = 0.0;
	for (auto vi : v) s += vi*vi;
	return sqrt(s);
}

double l2_sqrnorm(const std::vector<double>& v)
{
	double s = 0.0;
	for (auto vi : v) s += vi*vi;
	return s;
}

double l2_norm(double* x, int n)
{
	double s = 0.0;
	for (int i = 0; i < n; ++i) s += x[i]*x[i];
	return sqrt(s);
}
