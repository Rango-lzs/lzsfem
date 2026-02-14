#include "JFNKMatrix.h"
#include "femcore/NewtonSolver/NewtonSolver.h"
#include "femcore/FEModel.h"
#include "femcore/Domain/RgDomain.h"
#include "../FEMesh.h"

JFNKMatrix::JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K) : m_pns(pns), m_K(K)
{
	m_nrow = m_ncol = pns->m_neq;
	m_nsize = 0;

	m_bauto_eps = false;
	m_eps = 1e-6;

	m_policy = ZERO_PRESCRIBED_DOFS;

	// TODO: For contact problems we'll need some mechanism to change the array size
	m_v.resize(m_nrow);
	m_R.resize(m_nrow);

	// figure out the free and prescribed equation numbers
	m_freeDofs.clear();
	m_prescribedDofs.clear();

	FEModel* fem = m_pns->GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < node.m_dofs.size(); ++j)
		{
			int id = node.m_dofs[j];

			if (id >= 0) m_freeDofs.push_back(id);
			if (id < -1) m_prescribedDofs.push_back(-id - 2);
		}
	}

	// Add element dofs
	//Rango TODO
	/*for (int i = 0; i < mesh.Domains(); ++i)
	{
		RgDomain& dom = mesh.Domain(i);
		int NEL = dom.Elements();
		for (int j = 0; j < NEL; ++j)
		{
			RgElement& elj = dom.ElementRef(j);
			if (elj.m_lm >= 0) m_freeDofs.push_back(elj.m_lm);
		}
	}*/

	// make sure it all matches
	assert(m_freeDofs.size() + m_prescribedDofs.size() == m_pns->m_neq);
}

//! set Matrix policy
void JFNKMatrix::SetPolicy(MultiplyPolicy p)
{
	m_policy = p;
}

//! Set the forward difference epsilon
void JFNKMatrix::SetEpsilon(double eps)
{
	m_eps = eps;
}

//! Create a sparse Matrix from a sparse-Matrix profile
void JFNKMatrix::Create(SparseMatrixProfile& MP) 
{ 
	m_K->Create(MP); 
	m_nrow = m_K->Rows();
	m_ncol = m_K->Columns();
	m_nsize = m_K->NonZeroes(); 
}

void JFNKMatrix::SetReferenceResidual(std::vector<double>& R0)
{
	m_R0 = R0;
}

bool JFNKMatrix::mult_vector(double* x, double* r)
{
	int neq = 0;// (int)m_pns->m_ui.size();

	if (m_policy == ZERO_PRESCRIBED_DOFS)
	{
		for (int i = 0; i < m_freeDofs.size(); ++i)
		{
			int id = m_freeDofs[i];
			m_v[id] = x[id];
		}
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			m_v[id] = 0.0;
		}
	}
	else
	{
		for (int i = 0; i < m_freeDofs.size(); ++i)
		{
			int id = m_freeDofs[i];
			m_v[id] = 0.0;
		}
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			m_v[id] = x[id];
		}
	}

	double eps = m_eps;
	if (m_bauto_eps)
	{
		std::vector<double> u;
		m_pns->GetSolutionVector(u);

		assert(neq == Rows());
		double norm_v = l2_norm(m_v);

		eps = 0.0;
		if (norm_v != 0.0)
		{
			for (int i = 0; i < neq; ++i) eps += fabs(u[i]);
			eps *= m_eps / (neq*norm_v);
		}

		eps += m_eps;
	}

	// multiply by eps
	for (int i = 0; i < neq; ++i) m_v[i] *= eps;

	m_pns->Update2(m_v);
	if (m_pns->Residual(m_R) == false) return false;

	for (int i = 0; i < m_freeDofs.size(); ++i)
	{
		int id = m_freeDofs[i];
		r[id] = (m_R0[id] - m_R[id]) / eps;
	}

	if (m_policy == ZERO_PRESCRIBED_DOFS)
	{
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			r[id] = x[id];
		}
	}
	else
	{
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			r[id] = 0.0;
		}
	}

	return  true;
}
