#include "FEResidualVector.h"
#include "FERigidBody.h"
#include "FECore/DOFS.h"
#include "FEMechModel.h"
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEAnalysis.h>
#include "FESolidSolver2.h"
using namespace std;

//-----------------------------------------------------------------------------
FEResidualVector::FEResidualVector(FEModel& fem, vector<double>& R, vector<double>& Fr) : FEGlobalVector(fem, R, Fr)
{
}

//-----------------------------------------------------------------------------
FEResidualVector::~FEResidualVector()
{
}

//-----------------------------------------------------------------------------
void FEResidualVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe, bool bdom)
{
    
    vector<double>& R = m_R;
    
    int i, I, n;
    
    Vector3d a, d;
    
    {
        // assemble the element residual into the global residual
        int ndof = (int)fe.size();
        for (i=0; i<ndof; ++i)
        {
            
            I = elm[i];
            
            if ( I >= 0){
#pragma omp atomic
                R[I] += fe[i];
            }
            // TODO: Find another way to store reaction forces
            
            else if (-I-2 >= 0){
#pragma omp atomic
                m_Fr[-I-2] -= fe[i];
            }
        }
        
        
        int ndn = ndof / (int)en.size();
        // if there are linear constraints we need to apply them
        
        
		// process linear constraints
		FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
		if (LCM.LinearConstraints())
		{
			LCM.AssembleResidual(R, en, elm, fe);
		}
        
        // If there are rigid bodies we need to look for rigid dofs
		FEMechModel* fem = dynamic_cast<FEMechModel*>(&m_fem);
        if (fem && (fem->RigidBodies() > 0))
        {
            int *lm;
            for (i=0; i<ndof; i+=ndn)
            {
				int nid = en[i / ndn];
				if (nid >= 0)
				{
					FENode& node = m_fem.GetMesh().Node(nid);
					if (node.m_rid >= 0)
					{

						{
							Vector3d F(fe[i], fe[i + 1], fe[i + 2]);

							// this is an interface dof
							// get the rigid body this node is connected to
							FERigidBody& RB = *fem->GetRigidBody(node.m_rid);
							lm = RB.m_LM;

							// add to total torque of this body
							a = node.m_rt - RB.m_rt;
							Vector3d m = a ^ F;
							Vector3d f = F;

							// TODO: This code is only relevant when called from the shell domain residual and applies
							//	     the reaction of the back-face nodes.
							if (bdom)
							{
								if (node.HasFlags(FENode::SHELL) && node.HasFlags(FENode::RIGID_CLAMP)) {
									Vector3d d = node.m_dt;
									Vector3d b = a - d;
									Vector3d Fd(fe[i + 3], fe[i + 4], fe[i + 5]);
									f += Fd;
									m += b ^ Fd;
								}
							}

							n = lm[3];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += m.x;
							}
#pragma omp atomic
							RB.m_Mr.x -= m.x;
							n = lm[4];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += m.y;
							}

#pragma omp atomic
							RB.m_Mr.y -= m.y;
							n = lm[5];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += m.z;
							}
#pragma omp atomic
							RB.m_Mr.z -= m.z;
							/*
							 // if the rotational degrees of freedom are constrained for a rigid node
							 // then we need to add an additional component to the residual
							 if (node.m_ID[m_dofRU] == lm[3])
							 {
							 d = node.m_Dt;
							 n = lm[3]; if (n >= 0) R[n] += d.y*F.z-d.z*F.y; RB.m_Mr.x -= d.y*F.z-d.z*F.y;
							 n = lm[4]; if (n >= 0) R[n] += d.z*F.x-d.x*F.z; RB.m_Mr.y -= d.z*F.x-d.x*F.z;
							 n = lm[5]; if (n >= 0) R[n] += d.x*F.y-d.y*F.x; RB.m_Mr.z -= d.x*F.y-d.y*F.x;
							 }
							 */
							 // add to global force vector
							n = lm[0];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += f.x;
							}
#pragma omp atomic
							RB.m_Fr.x -= f.x;
							n = lm[1];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += f.y;
							}
#pragma omp atomic
							RB.m_Fr.y -= f.y;

							n = lm[2];
							if (n >= 0)
							{
#pragma omp atomic
								R[n] += f.z;
							}
#pragma omp atomic
							RB.m_Fr.z -= f.z;
						}
					}
				}
            }
        }
    }
}

//! Assemble into this global vector
void FEResidualVector::Assemble(int node_id, int dof, double f)
{
	// get the mesh
	FEMechModel& mechModel = dynamic_cast<FEMechModel&>(m_fem);
	FEMesh& mesh = m_fem.GetMesh();

	// get the equation number
	FENode& node = mesh.Node(node_id);
	int n = node.m_ID[dof];

	// assemble into global vector
	if (n >= 0) {
#pragma omp atomic
		m_R[n] += f;
	}
	else {
		FESolidSolver2* solver = dynamic_cast<FESolidSolver2*>(m_fem.GetCurrentStep()->GetFESolver());
		if (solver)
		{
			FERigidSolver* rigidSolver = solver->GetRigidSolver();
			rigidSolver->AssembleResidual(node_id, dof, f, m_R);
		}
	}
}
