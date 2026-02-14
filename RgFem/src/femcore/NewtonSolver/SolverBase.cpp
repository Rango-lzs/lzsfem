/*********************************************************************
 * \file   FESolver.cpp
 * \brief  Implementation of FESolver base class
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "femcore/NewtonSolver/SolverBase.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/Matrix/FEGlobalMatrix.h"
#include "femcore/FENode.h"
#include "logger/log.h"

DEFINE_META_CLASS(FESolver, FEObjectBase,"");

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FESolver, FEObjectBase)
    ADD_PARAMETER(m_bwopt, "optimize_bw");
    ADD_PARAMETER(m_msymm, "symmetric_stiffness");
    ADD_PARAMETER(m_neq, "equation_count");
END_PARAM_DEFINE();


//-----------------------------------------------------------------------------
// Scheme for assigning equation numbers
// STAGGERED: | a0, b0, a1, b1, ..., an, bn |
// BLOCK    : | a0, a1, ..., an, b0, b1, ..., bn |
enum EQUATION_SCHEME
{
    STAGGERED,
    BLOCK
};

//-----------------------------------------------------------------------------
FESolver::FESolver()
{
    m_bwopt = true;
    m_msymm = REAL_SYMMETRIC;
    m_neq = 0;
    
    // Statistics
    m_nrhs = 0;
    m_niter = 0;
    m_nref = 0;
    m_ntotref = 0;

    m_eq_scheme = 0;
}

//-----------------------------------------------------------------------------
FESolver::~FESolver()
{
}

//-----------------------------------------------------------------------------
void FESolver::Serialize(DumpStream& ar)
{
    FEObjectBase::Serialize(ar);
    
    ar & m_neq;
    ar & m_nrhs;
    ar & m_niter;
    ar & m_nref;
    ar & m_ntotref;
    ar & m_part;
    ar & m_dofMap;
}

//-----------------------------------------------------------------------------
bool FESolver::Init()
{
    FEModel* fem = GetFEModel();
    if (!fem) return false;
    
    // Initialize equation system
    if (!InitEquations()) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
void FESolver::Clean()
{
    // Clean up solver data
    m_Var.clear();
    m_dofMap.clear();
    m_part.clear();
}

//-----------------------------------------------------------------------------
void FESolver::Rewind()
{
    // Default implementation - do nothing
}

//-----------------------------------------------------------------------------
void FESolver::Reset()
{
    // Reset counters
    m_nrhs = 0;
    m_niter = 0;
    m_nref = 0;
    m_ntotref = 0;
}

//-----------------------------------------------------------------------------
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
        /*FENodeReorder mod;
        mod.Apply(mesh, P);*/
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
    if (m_eq_scheme == EQUATION_SCHEME::STAGGERED)  //[x0 y0  x1 y1 x2 y2]  x,y均为变量，自身可有多个分量
    {
        for (int i = 0; i < mesh.Nodes(); ++i)      // 先遍历节点
        {
            FENode& node = mesh.Node(P[i]);
            if (!node.HasFlags(FENode::EXCLUDE))
            {
                for (int jDof = 0; jDof < node.dofSize(); ++jDof)
                {
                    if (node.is_active(jDof))
                    {
                        int state = node.getDofState(jDof);
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
                        node.setDofIdx(jDof, -1);
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
        // DOFS& dofs = fem.GetDOFS();
        const RgDofSchema& dofs = fem.GetDofSchema();
        for (int jDof = 0; jDof < dofs.GetDofsPerNode(); ++jDof)
        {
            int neq0 = neq;
            for (int i = 0; i < mesh.Nodes(); ++i)
            {
                FENode& node = mesh.Node(P[i]);
                if (node.HasFlags(FENode::EXCLUDE) == false)
                {
                    if (node.is_active(jDof))
                    {
                        int bcl = node.getDofState(jDof);
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
            // assign partitions
            if (neq - neq0 > 0)
                m_part.push_back(neq - neq0);
        }
    }

    // store the number of equations
    m_neq = neq;
    assert(m_dofMap.size() == m_neq);
  
    // Build DOF map, what this purpose?
    m_dofMap.resize(m_neq);
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        int ndof = node.dofSize();
        
        for (int j = 0; j < ndof; ++j)
        {
            int eq = node.getDofIdx(j);
            if (eq >= 0)
            {
                m_dofMap[eq] = j;
            }
        }
    }
    
    RgLog("Equations allocated: %d\n", m_neq);  
    return true;
}

//-----------------------------------------------------------------------------
void FESolver::AddEquations(int neq, int partition)
{
    m_neq += neq;
}

//-----------------------------------------------------------------------------
void FESolver::SetPartitions(const std::vector<int>& part)
{
    m_part = part;
}

//-----------------------------------------------------------------------------
int FESolver::GetPartitionSize(int partition)
{
    if (m_part.empty() || partition < 0 || partition >= (int)m_part.size())
        return 0;
    
    return m_part[partition];
}

//-----------------------------------------------------------------------------
void FESolver::Update(std::vector<double>& u)
{
    // Default implementation - derived classes should override
}

//-----------------------------------------------------------------------------
int FESolver::MatrixSymmetryFlag() const
{
    return m_msymm;
}

//-----------------------------------------------------------------------------
MatrixType FESolver::MatType() const
{
    switch (m_msymm)
    {
        case REAL_SYMMETRIC: return MatrixType::REAL_SYMMETRIC;
        case REAL_UNSYMMETRIC: return MatrixType::REAL_UNSYMMETRIC;
        default: return MatrixType::REAL_UNSYMMETRIC;
    }
}

//-----------------------------------------------------------------------------
void FESolver::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
    FEModel* fem = GetFEModel();
    if (!fem) return;
    
    FEMesh& mesh = fem->GetMesh();
    
    // Build the matrix profile from mesh connectivity
    if (breset)
    {
        G.Clear();
        G.Create(mesh,m_neq, 0);
    }
    
    // Add entries for each domain
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.isActive())
        {
            dom.BuildMatrixProfile(G);
        }
    }
    
    // Add entries for contact
    for (int i = 0; i < fem->SurfacePairConstraints(); ++i)
    {
        /* FESurfacePairConstraint* pci = fem->SurfacePairConstraint(i);
         if (pci && pci->IsActive())
         {
             pci->BuildMatrixProfile(G);
         }*/
    }
    
    // Add entries for nonlinear constraints
    for (int i = 0; i < fem->NonlinearConstraints(); ++i)
    {
        /*FENLConstraint* plc = fem->NonlinearConstraint(i);
        if (plc && plc->IsActive())
        {
            plc->BuildMatrixProfile(G);
        }*/
    }
}

//-----------------------------------------------------------------------------
bool FESolver::HasActiveDofs(const FEDofList& dofs)
{
    for (const auto& var : m_Var)
    {
        if (var.m_dofs)
        {
            for (int i = 0; i < dofs.Size(); ++i)
            {
                if (var.m_dofs->Contains(dofs[i]))
                    return true;
            }
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
int FESolver::GetActiveDofMap(std::vector<int>& dofMap)
{
    dofMap = m_dofMap;
    return (int)dofMap.size();
}

//-----------------------------------------------------------------------------
FENodalDofInfo FESolver::GetDOFInfoFromEquation(int ieq)
{
    FENodalDofInfo info;
    
    if (ieq < 0 || ieq >= m_neq)
        return info;
    
    FEModel* fem = GetFEModel();
    if (!fem) return info;
    
    FEMesh& mesh = fem->GetMesh();
    DOFS& dofs = fem->GetDOFS();
    
    // Find the node and dof
    int NN = mesh.Nodes();
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        int ndof = node.dofSize();
        
        for (int j = 0; j < ndof; ++j)
        {
            if (node.getDofIdx(j) == ieq)
            {
                info.m_eq = ieq;
                info.m_node = i;
                info.m_dof = j;
                info.szdof = dofs.GetDOFName(j);
                return info;
            }
        }
    }
    
    return info;
}

//-----------------------------------------------------------------------------
double FESolver::ExtractSolutionNorm(const std::vector<double>& v, const FEDofList& dofs) const
{
    if (v.empty() || dofs.IsEmpty())
        return 0.0;
    
    FEModel* fem = GetFEModel();
    if (!fem) return 0.0;
    
    FEMesh& mesh = fem->GetMesh();
    
    double norm = 0.0;
    int NN = mesh.Nodes();
    
    for (int i = 0; i < NN; ++i)
    {
        const FENode& node = mesh.Node(i);
        
        for (int j = 0; j < dofs.Size(); ++j)
        {
            int dofId = dofs[j];
            if (dofId >= 0)
            {
                int eq = node.getDofIdx(dofId);
                if (eq >= 0 && eq < (int)v.size())
                {
                    norm += v[eq] * v[eq];
                }
            }
        }
    }
    
    return sqrt(norm);
}

//-----------------------------------------------------------------------------
std::vector<double> FESolver::GetSolutionVector() const
{
    std::vector<double> U(m_neq, 0.0);
    
    FEModel* fem = GetFEModel();
    if (!fem) return U;
    
    FEMesh& mesh = fem->GetMesh();
    int NN = mesh.Nodes();
    
    for (int i = 0; i < NN; ++i)
    {
        const FENode& node = mesh.Node(i);
        int ndof = node.dofSize();
        
        for (int j = 0; j < ndof; ++j)
        {
            int eq = node.getDofIdx(j);
            if (eq >= 0)
            {
                U[eq] = node.get(j);
            }
        }
    }
    
    return U;
}

//-----------------------------------------------------------------------------
void FESolver::AddSolutionVariable(FEDofList* dofs, int order, const char* szname)
{
    m_Var.push_back(FESolutionVariable(szname, dofs, order));
}
