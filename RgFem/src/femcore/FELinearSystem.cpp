#include "FELinearSystem.h"
#include "FELinearConstraintManager.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FELinearSystem::FELinearSystem(FESolver* solver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm) : m_K(K), m_F(F), m_u(u), m_solver(solver)
{
	m_bsymm = bsymm;
}

//-----------------------------------------------------------------------------
FELinearSystem::~FELinearSystem()
{

}

//-----------------------------------------------------------------------------
// get symmetry flag
bool FELinearSystem::IsSymmetric() const
{
	return m_bsymm;
}

//-----------------------------------------------------------------------------
// Get the solver that is using this linear system
FESolver* FELinearSystem::GetSolver()
{
	return m_solver;
}

//-----------------------------------------------------------------------------
//! assemble global stiffness Matrix
void FELinearSystem::Assemble(const FEElementMatrix& ke)
{
	if ((ke.rows() == 0) || (ke.columns() == 0)) return;

	// assemble into the global stiffness
	m_K.Assemble(ke);

	// check the prescribed contributions
	SparseMatrix& K = m_K;
	int N = ke.rows();
	int neq = m_K.Rows();

	// loop over columns
	const std::vector<int>& lmi = ke.RowIndices();
	const std::vector<int>& lmj = ke.ColumnsIndices();
	for (int j = 0; j<N; ++j)
	{
		int J = -lmj[j] - 2;
		if ((J >= 0) && (J<neq))
		{
			// dof j is a prescribed degree of freedom

			// loop over rows
			for (int i = 0; i<N; ++i)
			{
				int I = lmi[i];
				if (I >= 0)
				{
					// dof i is not a prescribed degree of freedom
#pragma omp atomic
					m_F[I] -= ke[i][j] * m_u[J];
				}
			}

			// set the diagonal element of K to 1
			K.set(J, J, 1);
		}
	}

#pragma omp critical
	{
	FEModel* fem = m_solver->GetFEModel();
	FELinearConstraintManager& LCM = fem->GetLinearConstraintManager();
	if (LCM.LinearConstraints())
	{
		const std::vector<int>& en = ke.Nodes();
		LCM.AssembleStiffness(m_K, m_F, m_u, en, lmi, lmj, ke);
	}
	} // omp critical
}

//-----------------------------------------------------------------------------
void FELinearSystem::AssembleRHS(std::vector<int>& lm, Matrix& ke, std::vector<double>& U)
{
	int ne = (int)lm.size();
	for (int j = 0; j<ne; ++j)
	{
		if (lm[j] >= 0)
		{
			double q = 0;
			for (int k = 0; k<ne; ++k)
			{
				if (lm[k] >= 0) q += ke[j][k] * U[lm[k]];
				else if (-lm[k] - 2 >= 0) q += ke[j][k] * U[-lm[k] - 2];
			}
			m_F[lm[j]] += q;
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearSystem::AssembleRHS(std::vector<int>& lm, std::vector<double>& fe)
{
	const int n = (int)lm.size();
	for (int i = 0; i<n; ++i)
	{
		int nid = lm[i];
		if (nid >= 0) 
		{
			m_F[nid] += fe[i];
		}
	}
}
