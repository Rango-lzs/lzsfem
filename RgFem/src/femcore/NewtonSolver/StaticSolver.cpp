/*********************************************************************
 * \file   FEStaticSolver.cpp
 * \brief  Implementation of static solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "femcore/NewtonSolver/StaticSolver.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FELinearSystem.h"
#include "femcore/Domain/RgDomain.h"
#include "logger/log.h"
#include "../BoundaryCondition/RgBoundaryCondition.h"

DEFINE_META_CLASS(FEStaticSolver, FENewtonSolver,"static");

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEStaticSolver, FENewtonSolver)
    ADD_PARAMETER(m_Dtol, "dtol");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEStaticSolver::FEStaticSolver()
{
    m_Dtol = 0.001;  // Default displacement tolerance
    
    // Add solution variable for displacements
    // This will be set up properly in Init()
}

//-----------------------------------------------------------------------------
FEStaticSolver::~FEStaticSolver()
{
}

//-----------------------------------------------------------------------------
void FEStaticSolver::Serialize(DumpStream& ar)
{
    FENewtonSolver::Serialize(ar);
    
    ar & m_Dtol;
    ar & m_Fr;
    ar & m_Fint;
    ar & m_Fext;
}

//-----------------------------------------------------------------------------
bool FEStaticSolver::Init()
{
    // Call base class initialization
    if (!FENewtonSolver::Init() == false)
        return false;

    FEModel& fem = *GetFEModel();
    
    // Get displacement DOFs
    DOFS& dofs = fem.GetDOFS();
    m_dofU.Clear();
    m_dofU.AddDof(dofs.GetDOF("x"));
    m_dofU.AddDof(dofs.GetDOF("y"));
    m_dofU.AddDof(dofs.GetDOF("z"));
    
    // Add displacement as a solution variable
    AddSolutionVariable(&m_dofU, 2, "displacement", m_Dtol);
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEStaticSolver::InitEquations()
{
    // Initialize equation system
    return FENewtonSolver::InitEquations();
}

//-----------------------------------------------------------------------------
bool FEStaticSolver::InitStep(double time)
{
    FEModel& fem = *GetFEModel();
    
    // Update time
    FETimeInfo& tp = fem.GetTime();
    tp.currentTime = time;
    
    // Allocate force vectors if needed
    int neq = m_neq;
    if (m_Fr.size() != neq)
    {
        m_Fr.assign(neq, 0.0);
        m_Fint.assign(neq, 0.0);
        m_Fext.assign(neq, 0.0);
    }
    
    // Zero force vectors
    std::fill(m_Fr.begin(), m_Fr.end(), 0.0);
    std::fill(m_Fint.begin(), m_Fint.end(), 0.0);
    std::fill(m_Fext.begin(), m_Fext.end(), 0.0);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEStaticSolver::PrepStep()
{
    // Zero the total displacements
    std::fill(m_Ui.begin(), m_Ui.end(), 0.0);
    
    // Call base class
    FENewtonSolver::PrepStep();
}

//-----------------------------------------------------------------------------
bool FEStaticSolver::StiffnessMatrix()
{
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    
    // Get the linear system
    FELinearSystem LS(this, *m_pK, m_Fd, m_ui, false);
    
    // Assemble stiffness from all domains
    FEMesh& mesh = fem.GetMesh();
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.isActive())
        {
            //dom.StiffnessMatrix(LS);
        }
    }
    
    // Add contact stiffness
    ContactStiffness(LS);
    
    // Add nonlinear constraint stiffness
    NonLinearConstraintStiffness(LS, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEStaticSolver::Residual(std::vector<double>& R)
{
    // Zero the residual
    std::fill(R.begin(), R.end(), 0.0);
    
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    
    // Create global vector
    FEGlobalVector RHS(fem, R, m_Fd);
    
    // Internal forces (negative because we solve K*u = F - Fint)
    InternalForces(RHS);
    
    // External forces
    ExternalForces(RHS);
    
    // Contact forces
    ContactForces(RHS);
    
    // Nonlinear constraint forces
    NonLinearConstraintForces(RHS, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEStaticSolver::InternalForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // Assemble internal forces from all domains
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.isActive())
        {
            //dom.InternalForces(R);
        }
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::ExternalForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
    // Apply boundary conditions
    int nbc = fem.BoundaryConditions();
    for (int i = 0; i < nbc; ++i)
    {
        RgBoundaryCondition* pbc = fem.BoundaryCondition(i);
        if (pbc && pbc->IsActive())
        {
            // pbc->LoadVector(R);
        }
    }
    
    // Apply model loads
    int nml = fem.ModelLoads();
    for (int i = 0; i < nml; ++i)
    {
        RgLoad* pml = fem.ModelLoad(i);
        //if (pml && pml->IsActive())
        //{
        //    // pml->LoadVector(R);
        //}
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::ContactForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
    // Contact interface forces
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        //FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
        //if (pci && pci->IsActive())
        //{
        //    // pci->LoadVector(R);
        //}
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    // Nonlinear constraint forces
    int nnlc = fem.NonlinearConstraints();
    for (int i = 0; i < nnlc; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        //if (plc && plc->IsActive())
        //{
        //    // plc->LoadVector(R, tp);
        //}
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::ContactStiffness(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    
    // Contact interface stiffness
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
        //if (pci && pci->IsActive())
        //{
        //    // pci->StiffnessMatrix(LS);
        //}
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    // Nonlinear constraint stiffness
    int nnlc = fem.NonlinearConstraints();
    for (int i = 0; i < nnlc; ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        //if (plc && plc->IsActive())
        //{
        //    // plc->StiffnessMatrix(LS, tp);
        //}
    }
}

//-----------------------------------------------------------------------------
void FEStaticSolver::Update(std::vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // Update nodal positions
    int N = mesh.Nodes();
    for (int i = 0; i < N; ++i)
    {
        FENode& node = mesh.Node(i);
        
        // Get displacement DOF indices
        int nx = node.getDofIdx(m_dofU[0]);
        int ny = node.getDofIdx(m_dofU[1]);
        int nz = node.getDofIdx(m_dofU[2]);
        
        // Update position
        if (nx >= 0) node.m_rt.x = node.m_r0.x + node.get(m_dofU[0]) + ui[nx];
        if (ny >= 0) node.m_rt.y = node.m_r0.y + node.get(m_dofU[1]) + ui[ny];
        if (nz >= 0) node.m_rt.z = node.m_r0.z + node.get(m_dofU[2]) + ui[nz];
        
        // Update displacement values
        if (nx >= 0) node.set(m_dofU[0], node.m_rt.x - node.m_r0.x);
        if (ny >= 0) node.set(m_dofU[1], node.m_rt.y - node.m_r0.y);
        if (nz >= 0) node.set(m_dofU[2], node.m_rt.z - node.m_r0.z);
    }
    
    // Update mesh
    mesh.Update(fem.GetTime());
}
