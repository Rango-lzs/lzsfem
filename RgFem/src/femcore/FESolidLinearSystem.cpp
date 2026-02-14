#include "FESolidLinearSystem.h"
#include "femcore/FELinearConstraintManager.h"
#include "femcore/FEModel.h"

FESolidLinearSystem::FESolidLinearSystem(FESolver* solver, FERigidSolver* rigidSolver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm, double alpha, int nreq) : FELinearSystem(solver, K, F, u, bsymm)
{
	m_rigidSolver = rigidSolver;
	m_alpha = alpha;
	m_nreq = nreq;
	m_stiffnessScale = 1.0;
}

// scale factor for stiffness Matrix
void FESolidLinearSystem::StiffnessAssemblyScaleFactor(double a)
{
	m_stiffnessScale = a;
}

void FESolidLinearSystem::Assemble(const FEElementMatrix& ke)
{
	// Rigid joints require a different assembly approach in that we can do 
	// a direct assembly as defined by the base class. 
	// Currently, we assume that if the node list of the element Matrix is not
	// defined, then we are dealing with rigid joints.
	if (ke.Nodes().empty())
	{
		FELinearSystem::Assemble(ke);
	}
	else
	{
		// assemble into global stiffness Matrix
		if (m_stiffnessScale == 1.0)
		{
			m_K.Assemble(ke);
		}
		else
		{
			// NOTE: What is doing here?
			FEElementMatrix kes(ke, m_stiffnessScale);
			m_K.Assemble(kes);
		}

		// get the std::vector that stores the prescribed BC values
		std::vector<double>& ui = m_u;

		// adjust for linear constraints
		FEModel* fem = m_solver->GetFEModel();
		FELinearConstraintManager& LCM = fem->GetLinearConstraintManager();
		if (LCM.LinearConstraints() > 0)
		{
			#pragma omp critical 
			LCM.AssembleStiffness(m_K, m_F, m_u, ke.Nodes(), ke.RowIndices(), ke.ColumnsIndices(), ke);
		}

		// adjust stiffness Matrix for prescribed degrees of freedom
		// NOTE: I had to comment this if statement out since otherwise
		//       poroelastic DOF's that are set as free-draining in the
		//       sliding2 contact code are skipt and zeroes will appear
		//       on the diagonal of the stiffness Matrix.
		//	if (m_fem.m_DC.size() > 0)
		{
			SparseMatrix& K = m_K;

			int N = ke.rows();

			// loop over columns
			const std::vector<int>& elmi = ke.RowIndices();
			const std::vector<int>& elmj = ke.ColumnsIndices();
			for (int j = 0; j < N; ++j)
			{
				int J = -elmj[j] - 2;
				if ((J >= 0) && (J < m_nreq))
				{
					// dof j is a prescribed degree of freedom

					// loop over rows
					for (int i = 0; i < N; ++i)
					{
						int I = elmi[i];
						if (I >= 0)
						{
							// dof i is not a prescribed degree of freedom
							#pragma omp atomic
							m_F[I] -= ke[i][j] * ui[J];
						}
					}

					// set the diagonal element of K to 1
					K.set(J, J, 1);
				}
			}
		}

		// see if there are any rigid body dofs here
		#pragma omp critical 
		/*m_rigidSolver->RigidStiffness(m_K, m_u, m_F, ke, m_alpha);*/
	}
}
