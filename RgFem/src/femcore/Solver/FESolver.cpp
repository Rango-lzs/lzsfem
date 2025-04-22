/*****************************************************************
 * \file   FESolver.cpp
 * \brief
 *
 * \author 11914
 * \date   March 2025
 *********************************************************************/

#include "FESolver.h"

#include "basicio/DumpStream.h"
#include "femcore/FELinearConstraintManager.h"
#include "femcore/Domain/FEDomain.h"
#include "femcore/FEMesh.h"
#include "femcore/FEModel.h"
#include "femcore/FENode.h"
#include "femcore/FENLConstraint.h"
#include "femcore/FENodalLoad.h"
#include "FENodeReorder.h"
#include "femcore/FESurfacePairConstraint.h"
#include "LinearSolver.h"

BEGIN_PARAM_DEFINE(FESolver, FEObjectBase)
BEGIN_PARAM_GROUP("linear system");
ADD_PARAMETER(m_msymm, "symmetric_stiffness", 0, "non-symmetric\0symmetric\0symmetric structure\0");
ADD_PARAMETER(m_eq_scheme, "equation_scheme", 0, "staggered\0block\0");
ADD_PARAMETER(m_eq_order, "equation_order", 0, "default\0reverse\0febio2\0");
ADD_PARAMETER(m_bwopt, "optimize_bw");
END_PARAM_GROUP();
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FESolver::FESolver(FEModel* fem)
    : FEObjectBase(fem)
{
    m_niter = 0;
    m_nref = 0;
    m_baugment = false;
    m_naug = 0;
    m_neq = 0;
    m_bwopt = false;
}

//-----------------------------------------------------------------------------
FESolver::~FESolver()
{
}

//-----------------------------------------------------------------------------
void FESolver::SetEquationScheme(int scheme)
{
    m_eq_scheme = scheme;
}

//-----------------------------------------------------------------------------
//! set the linear system partitions
void FESolver::SetPartitions(const std::vector<int>& part)
{
    m_part = part;
}

//-----------------------------------------------------------------------------
//! Get the size of a partition
int FESolver::GetPartitionSize(int partition)
{
    assert((partition >= 0) && (partition < (int)m_part.size()));
    if ((partition >= 0) && (partition < (int)m_part.size()))
        return m_part[partition];
    else
        return 0;
}

//-----------------------------------------------------------------------------
void FESolver::Clean()
{
}

//-----------------------------------------------------------------------------
void FESolver::Reset()
{
    m_niter = 0;
    m_nref = 0;
    m_naug = 0;
}

//-----------------------------------------------------------------------------
// get the linear solver
LinearSolver* FESolver::GetLinearSolver()
{
    return nullptr;
}

//-----------------------------------------------------------------------------
//! Matrix symmetry flag
int FESolver::MatrixSymmetryFlag() const
{
    return m_msymm;
}

//-----------------------------------------------------------------------------
//! get matrix type
MatrixType FESolver::MatType() const
{
    return static_cast<MatrixType>(m_msymm);
}

//-----------------------------------------------------------------------------
// extract the (square) norm of a solution std::vector
double FESolver::ExtractSolutionNorm(const std::vector<double>& v, const FEDofList& dofs) const
{
    assert(v.size() == m_dofMap.size());
    double norm = 0;
    for (int n = 0; n < dofs.Size(); ++n)
    {
        for (int i = 0; i < v.size(); ++i)
        {
            if (m_dofMap[i] == dofs[n])
                norm += v[i] * v[i];
        }
    }
    return norm;
}

//-----------------------------------------------------------------------------
// see if the dofs in the dof list are active in this solver
bool FESolver::HasActiveDofs(const FEDofList& dof)
{
    assert(dof.IsEmpty() == false);
    if (dof.IsEmpty())
        return true;

    assert(m_Var.size());
    for (int i = 0; i < dof.Size(); ++i)
    {
        int dof_i = dof[i];

        for (int i = 0; i < m_Var.size(); ++i)
        {
            FESolutionVariable& vi = m_Var[i];
            if (vi.m_dofs->Contains(dof_i))
            {
                return true;
            }
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
// get the active dof map (returns nr of functions)
int FESolver::GetActiveDofMap(std::vector<int>& activeDofMap)
{
    // get the dof map
    int neq = (int)m_dofMap.size();
    if (m_dofMap.empty() || (m_dofMap.size() < neq))
        return -1;

    // We need the partitions here, but for now we assume that
    // it is the first partition

    // The dof map indices point to the dofs as defined by the variables.
    // Since there could be more dofs than actually used in the linear system
    // we need to reindex this map.
    // First, find the min and max
    int imin = m_dofMap[0], imax = m_dofMap[0];
    for (size_t i = 0; i < neq; ++i)
    {
        if (m_dofMap[i] > imax)
            imax = m_dofMap[i];
        if (m_dofMap[i] < imin)
            imin = m_dofMap[i];
    }

    // create the conversion table
    int nsize = imax - imin + 1;
    std::vector<int> LUT(nsize, -1);
    for (size_t i = 0; i < neq; ++i)
    {
        LUT[m_dofMap[i] - imin] = 1;
    }

    // count how many dofs are actually used
    int nfunc = 0;
    for (size_t i = 0; i < nsize; ++i)
    {
        if (LUT[i] != -1)
            LUT[i] = nfunc++;
    }

    // now, reindex the dof map
    // allocate dof map
    activeDofMap.resize(neq);
    for (size_t i = 0; i < neq; ++i)
    {
        activeDofMap[i] = LUT[m_dofMap[i] - imin];
    }

    return nfunc;
}

//-----------------------------------------------------------------------------
//! build the matrix profile
void FESolver::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    DOFS& fedofs = fem.GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();

    // when reset is true we build the entire matrix profile
    // (otherwise we only build the "dynamic" profile)
    if (breset)
    {
        std::vector<int> elm;

        // Add all elements to the profile
        // Loop over all active domains
        for (int nd = 0; nd < mesh.Domains(); ++nd)
        {
            FEDomain& d = mesh.Domain(nd);
            d.BuildMatrixProfile(G);
        }

        // linear constraints
        FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
        LCM.BuildMatrixProfile(G);
    }
    else
    {
        // Do the "dynamic" profile. That is the part of the profile that always changes
        // This is mostly contact
        // do the nonlinear constraints
        int M = fem.NonlinearConstraints();
        for (int m = 0; m < M; ++m)
        {
            FENLConstraint* pnlc = fem.NonlinearConstraint(m);
            if (pnlc->IsActive())
                pnlc->BuildMatrixProfile(G);
        }

        // All following "elements" are nonstatic. That is, they can change
        // connectivity between calls to this function. All of these elements
        // are related to contact analysis (at this point).
        if (fem.SurfacePairConstraints() > 0)
        {
            // Add all contact interface elements
            for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
            {
                FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
                if (pci->IsActive())
                    pci->BuildMatrixProfile(G);
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! This function is called right before SolveStep and should be used to initialize
//! time dependent information and other settings.
bool FESolver::InitStep(double time)
{
    FEModel& fem = *GetFEModel();

    //// evaluate load controllers values at current time
    // fem.EvaluateLoadControllers(time);
    //// evaluate data generators at current time
    // fem.EvaluateDataGenerators(time);
    //// evaluate load parameters
    // fem.EvaluateLoadParameters();

    // re-validate materials
    // This is necessary since the material parameters can have changed (e.g. via load curves) and thus
    // a new validation needs to be done to see if the material parameters are still valid.
    if (fem.ValidateMaterials() == false)
        return false;

    return true;
}

// 方程编号
bool FESolver::InitEquations()
{
    // get the mesh
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    // clear partitions
    m_part.clear();

    // reorder the node numbers
    int NN = mesh.Nodes();
    std::vector<int> P(NN);

    // see if we need to optimize the bandwidth
    if (m_bwopt)
    {
        FENodeReorder mod;
        mod.Apply(mesh, P);
    }
    else
        for (int i = 0; i < NN; ++i)
            P[i] = i;

    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(P[i]);
        if (node.HasFlags(FENode::EXCLUDE))
            for (int j = 0; j < node.dofSize(); ++j)
                node.setDofIdx(j, -1);
    }
    m_dofMap.clear();

    // assign equations based on allocation scheme
    int neq = 0;
    if (m_eq_scheme == EQUATION_SCHEME::STAGGERED) //[x0 y0  x1 y1 x2 y2]  x,y均为变量，自身可有多个分量
    {
        DOFS& dofs = fem.GetDOFS();
        for (int i = 0; i < mesh.Nodes(); ++i) //先遍历节点
        {
            FENode& node = mesh.Node(P[i]);
            if (!node.HasFlags(FENode::EXCLUDE))
            {
                for (int iVar = 0; iVar < dofs.Variables(); ++iVar) //变量个数
                {
                    int nVar = dofs.GetVariableSize(iVar); //变量分量
                    for (int j = 0; j < nVar; ++j)
                    {
                        int jDof = dofs.GetDOF(iVar, j);
                        if (node.is_active(jDof))
                        {
                            int state = node.get_bc(jDof);
                            if (state == DOF_OPEN)
                            {
                                node.setDofIdx(jDof, neq++);
                                m_dofMap.push_back(jDof);
                            }
                            else if (state == DOF_FIXED)
                            {
                                node.setDofIdx(jDof, -1);
                            }
                            else if (state == DOF_PRESCRIBED)
                            {
                                node.setDofIdx(jDof, -neq - 2);
                                neq++;
                                m_dofMap.push_back(jDof);
                            }
                            else
                            {
                                assert(false);
                                return false;
                            }
                        }
                        else
                            node.setDofIdx(jDof,-1);
                    }
                }
            }
        }

        // assign partition
        m_part.push_back(neq);
    }
    else
    {
        // Assign equations numbers in blocks [x0 x1 x2 ……  y0 y1 y2]
        assert(m_eq_scheme == EQUATION_SCHEME::BLOCK);
        DOFS& dofs = fem.GetDOFS();
        for (int iVar = 0; iVar < dofs.Variables(); ++iVar)
        {
            int neq0 = neq;
            for (int i = 0; i < mesh.Nodes(); ++i)
            {
                FENode& node = mesh.Node(P[i]);
                if (node.HasFlags(FENode::EXCLUDE) == false)
                {
                    int nVar = dofs.GetVariableSize(iVar);
                    for (int j = 0; j < nVar; ++j)
                    {
                        int jDof = dofs.GetDOF(iVar, j);

                        if (node.is_active(jDof))
                        {
                            int bcl = node.get_bc(jDof);
                            if (bcl == DOF_FIXED)
                            {
                                node.setDofIdx(jDof, -1);
                            }
                            else if (bcl == DOF_OPEN)
                            {
                                node.setDofIdx(jDof, neq++);
                                m_dofMap.push_back(jDof);
                            }
                            else if (bcl == DOF_PRESCRIBED)
                            {
                                node.setDofIdx(jDof, -neq - 2);
                                neq++;
                                m_dofMap.push_back(jDof);
                            }
                            else
                            {
                                assert(false);
                                return false;
                            }
                        }
                        else
                            node.setDofIdx(jDof, -1);
                    }
                }
            }

            // assign partitions
            if (neq - neq0 > 0)
                m_part.push_back(neq - neq0);
        }
    }

    // store the number of equations
    m_neq = neq;

    assert(m_dofMap.size() == m_neq);

    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
//! add equations
void FESolver::AddEquations(int neq, int partition)
{
    m_neq += neq;
    m_part[partition] += neq;
}

//-----------------------------------------------------------------------------
void FESolver::Serialize(DumpStream& ar)
{
    FEObjectBase::Serialize(ar);
    ar& m_nrhs& m_niter& m_nref& m_ntotref& m_naug;
}

//-----------------------------------------------------------------------------
//! Update the state of the model
void FESolver::Update(std::vector<double>& u)
{
    assert(false);
};

//-----------------------------------------------------------------------------
// The augmentation is done after a time step converges and gives model components
// an opportunity to modify the model's state. This will usually require that the time
// step is solved again.
bool FESolver::Augment()
{
    FEModel& fem = *GetFEModel();

    const FETimeInfo& tp = fem.GetTime();

    // Assume we will pass (can't hurt to be optimistic)
    bool bconv = true;

    // Do contact augmentations
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
        if (pci->IsActive())
            bconv = (pci->Augment(m_naug, tp) && bconv);
    }

    // do nonlinear constraint augmentations
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc->IsActive())
            bconv = plc->Augment(m_naug, tp) && bconv;
    }

    // do domain augmentations
    FEMesh& mesh = fem.GetMesh();
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        bconv = dom.Augment(m_naug) && bconv;
    }

    fem.GetTime().augmentation++;

    return bconv;
}

//-----------------------------------------------------------------------------
// return the node (mesh index) from an equation number
FENodalDofInfo FESolver::GetDOFInfoFromEquation(int ieq)
{
    FENodalDofInfo info;
    info.m_eq = ieq;
    info.m_node = -1;
    info.m_dof = -1;
    info.szdof = "";

    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        const std::vector<int>& id = node.getDofs();
        for (int j = 0; j < id.size(); ++j)
        {
            if (id[j] == ieq)
            {
                info.m_node = node.GetID();
                info.m_dof = j;
                DOFS& Dofs = GetFEModel()->GetDOFS();
                info.szdof = Dofs.GetDOFName(info.m_dof);
                if (info.szdof == nullptr)
                    info.szdof = "???";
                return info;
            }
        }
    }
    return info;
}
