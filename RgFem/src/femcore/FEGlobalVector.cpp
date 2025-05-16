#include "FEGlobalVector.h"
#include "datastructure/Vector3d.h"
#include "FEModel.h"
#include "FENode.h"
#include "femcore/FEMesh.h"

//-----------------------------------------------------------------------------
FEGlobalVector::FEGlobalVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr) : m_fem(fem), m_R(R), m_Fr(Fr)
{
}

//-----------------------------------------------------------------------------
FEGlobalVector::~FEGlobalVector()
{

}

//-----------------------------------------------------------------------------
void FEGlobalVector::Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom)
{
	std::vector<double>& R = m_R;

	// assemble the element residual into the global residual
	int ndof = (int)fe.size();
	for (int i=0; i<ndof; ++i)
	{
		int I = elm[i];
		if ( I >= 0) {
#pragma omp atomic
			R[I] += fe[i];
		}
// TODO: Find another way to store reaction forces
		else if (-I-2 >= 0) {
#pragma omp atomic
			m_Fr[-I-2] -= fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! \todo This function does not add to m_Fr. Is this a problem?
void FEGlobalVector::Assemble(std::vector<int>& lm, std::vector<double>& fe)
{
	std::vector<double>& R = m_R;
	const int n = (int) lm.size();
	for (int i=0; i<n; ++i)
	{
		int nid = lm[i];
		if (nid >= 0) {
#pragma omp atomic
			R[nid] += fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! assemble a nodel value
void FEGlobalVector::Assemble(int nodeId, int dof, double f)
{
	// get the equation number
	FENode& node = m_fem.GetMesh().Node(nodeId);
	int n = node.m_dofs[dof];

	// assemble into global std::vector
	if (n >= 0) {
#pragma omp atomic
		m_R[n] += f;
	}
}
