/*********************************************************************
 * \file   FEImplicitDynamicSolver.cpp
 * \brief  Implementation of implicit dynamic solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "FEImplicitDynamicSolver.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FELinearSystem.h"
#include "femcore/Domain/RgDomain.h"
#include "logger/log.h"

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEImplicitDynamicSolver, FENewtonSolver)
    ADD_PARAMETER(m_Dtol, "dtol");
    ADD_PARAMETER(m_beta, "beta");
    ADD_PARAMETER(m_gamma, "gamma");
    ADD_PARAMETER(m_rhoi, "rhoi");
    ADD_PARAMETER(m_alphaf, "alphaf");
    ADD_PARAMETER(m_alpham, "alpham");
    ADD_PARAMETER(m_useDamping, "use_damping");
    ADD_PARAMETER(m_dampingRatio, "damping_ratio");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEImplicitDynamicSolver::FEImplicitDynamicSolver()
{
    // Default to Newmark average acceleration
    m_timeScheme = NEWMARK;
    m_beta = 0.25;
    m_gamma = 0.5;
    
    // Generalized-alpha defaults
    m_rhoi = 0.5;
    m_alphaf = 0.0;
    m_alpham = 0.0;
    
    // Convergence
    m_Dtol = 0.001;
    
    // Damping
    m_useDamping = false;
    m_dampingRatio = 0.0;
}

//-----------------------------------------------------------------------------
FEImplicitDynamicSolver::~FEImplicitDynamicSolver()
{
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::Serialize(DumpStream& ar)
{
    FENewtonSolver::Serialize(ar);
    
    ar & m_Dtol;
    ar & m_beta & m_gamma;
    ar & m_rhoi & m_alphaf & m_alpham;
    ar & m_Fr & m_Fint & m_Fext & m_Fin & m_Fd;
    ar & m_Uip;
}

//-----------------------------------------------------------------------------
bool FEImplicitDynamicSolver::Init()
{
    if (!FENewtonSolver::Init()) return false;
    
    FEModel& fem = *GetFEModel();
    DOFS& dofs = fem.GetDOFS();
    
    // Get DOFs
    m_dofU.Clear();
    m_dofU.AddDof(dofs.GetDOF("x"));
    m_dofU.AddDof(dofs.GetDOF("y"));
    m_dofU.AddDof(dofs.GetDOF("z"));
    
    m_dofV.Clear();
    m_dofV.AddDof(dofs.GetDOF("vx"));
    m_dofV.AddDof(dofs.GetDOF("vy"));
    m_dofV.AddDof(dofs.GetDOF("vz"));
    
    // Add solution variables
    AddSolutionVariable(&m_dofU, 2, "displacement", m_Dtol);
    AddSolutionVariable(&m_dofV, 2, "velocity", m_Dtol);
    
    // Initialize time integration
    InitializeTimeIntegration();
    
    // Allocate vectors
    int neq = m_neq;
    m_Fr.assign(neq, 0.0);
    m_Fint.assign(neq, 0.0);
    m_Fext.assign(neq, 0.0);
    m_Fin.assign(neq, 0.0);
    m_Fd.assign(neq, 0.0);
    m_Uip.assign(neq, 0.0);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::InitializeTimeIntegration()
{
    switch (m_timeScheme)
    {
        case NEWMARK:
            // Parameters already set
            RgLog("Using Newmark integration (beta=%.3f, gamma=%.3f)\n", m_beta, m_gamma);
            break;
            
        case GENERALIZED_ALPHA:
        {
            // Calculate alpha parameters from spectral radius
            m_alpham = (2.0 * m_rhoi - 1.0) / (m_rhoi + 1.0);
            m_alphaf = m_rhoi / (m_rhoi + 1.0);
            m_beta = 0.25 * (1.0 - m_alpham + m_alphaf) * (1.0 - m_alpham + m_alphaf);
            m_gamma = 0.5 - m_alpham + m_alphaf;
            
            RgLog("Using Generalized-Alpha (rhoi=%.3f, alphaf=%.3f, alpham=%.3f)\n", 
                  m_rhoi, m_alphaf, m_alpham);
        }
        break;
        
        case HHT_ALPHA:
        {
            // HHT-alpha method
            double alpha = m_alphaf;  // Use alphaf as HHT alpha
            m_alpham = 0.0;
            m_alphaf = alpha;
            m_beta = 0.25 * (1.0 + alpha) * (1.0 + alpha);
            m_gamma = 0.5 + alpha;
            
            RgLog("Using HHT-Alpha (alpha=%.3f)\n", alpha);
        }
        break;
    }
}

//-----------------------------------------------------------------------------
bool FEImplicitDynamicSolver::InitEquations()
{
    return FENewtonSolver::InitEquations();
}

//-----------------------------------------------------------------------------
bool FEImplicitDynamicSolver::InitStep(double time)
{
    FEModel& fem = *GetFEModel();
    
    // Update time
    FETimeInfo& tp = fem.GetTime();
    tp.currentTime = time;
    
    // Zero force vectors
    std::fill(m_Fr.begin(), m_Fr.end(), 0.0);
    std::fill(m_Fint.begin(), m_Fint.end(), 0.0);
    std::fill(m_Fext.begin(), m_Fext.end(), 0.0);
    std::fill(m_Fin.begin(), m_Fin.end(), 0.0);
    std::fill(m_Fd.begin(), m_Fd.end(), 0.0);
    
    // Predict values for new time step
    PredictValues();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::PredictValues()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    int NN = mesh.Nodes();
    
    // Newmark predictor
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofU = m_dofU[j];
            int dofV = m_dofV[j];
            
            // Get previous values
            double un = node.get(dofU);
            double vn = node.get(dofV);
            double an = node.get_prev(dofV) / dt;  // Approximate
            
            // Predict displacement
            double up = un + dt * vn + (0.5 - m_beta) * dt * dt * an;
            node.set(dofU, up);
            
            // Predict velocity
            double vp = vn + (1.0 - m_gamma) * dt * an;
            node.set(dofV, vp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::SolverWarnings()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Check time step size
    if (dt <= 0.0)
    {
        RgLogWarning("Time step size is zero or negative");
    }
    
    // Check for unstable parameters
    if (m_gamma < 0.5)
    {
        RgLogWarning("Gamma < 0.5 may lead to unconditionally unstable scheme");
    }
    
    if (m_beta < 0.25 * m_gamma * m_gamma)
    {
        RgLogWarning("Beta parameter may lead to unstable scheme");
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::PrepStep()
{
    FENewtonSolver::PrepStep();
    
    // Save previous displacement increment
    m_Uip = m_Ui;
}

//-----------------------------------------------------------------------------
bool FEImplicitDynamicSolver::StiffnessMatrix()
{
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Create linear system
    FELinearSystem LS(this, m_pK, m_Fd, m_ui, m_breshape);
    
    // Add mass matrix contribution: M / (beta * dt^2)
    double mass_coeff = 1.0 / (m_beta * dt * dt);
    MassMatrix(LS);
    LS.Scale(mass_coeff);
    
    // Add damping matrix contribution if enabled
    if (m_useDamping)
    {
        double damp_coeff = m_gamma / (m_beta * dt);
        DampingMatrix(LS);
        LS.Scale(damp_coeff);
    }
    
    // Add stiffness from domains
    FEMesh& mesh = fem.GetMesh();
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
        {
            dom.StiffnessMatrix(LS);
        }
    }
    
    // Add contact stiffness
    ContactStiffness(LS);
    
    // Add nonlinear constraint stiffness
    NonLinearConstraintStiffness(LS, tp);
    
    m_breshape = false;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEImplicitDynamicSolver::Residual(std::vector<double>& R)
{
    std::fill(R.begin(), R.end(), 0.0);
    
    FEModel& fem = *GetFEModel();
    const FETimeInfo& tp = fem.GetTime();
    
    FEGlobalVector RHS(fem, R, m_Fd);
    
    // Internal forces (negative)
    InternalForces(RHS);
    
    // External forces (positive)
    ExternalForces(RHS);
    
    // Inertial forces (negative)
    InertialForces(RHS);
    
    // Damping forces (negative)
    if (m_useDamping)
    {
        DampingForces(RHS);
    }
    
    // Contact forces
    ContactForces(RHS);
    
    // Nonlinear constraints
    NonLinearConstraintForces(RHS, tp);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::Update(std::vector<double>& ui)
{
    // Update kinematics
    UpdateKinematics(ui);
    
    // Update increments
    UpdateIncrements(m_Ui, ui);
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::UpdateKinematics(std::vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    int NN = mesh.Nodes();
    
    // Update using Newmark formulas
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofU = m_dofU[j];
            int dofV = m_dofV[j];
            
            int eq = node.m_ID[dofU];
            if (eq >= 0)
            {
                // Get current values
                double un = node.get(dofU);
                double vn = node.get(dofV);
                
                // Update displacement
                double du = ui[eq];
                double up = un + du;
                node.set(dofU, up);
                node.m_rt = node.m_r0 + Vector3d(
                    node.get(m_dofU[0]),
                    node.get(m_dofU[1]),
                    node.get(m_dofU[2])
                );
                
                // Update velocity using Newmark
                double dv = m_gamma / (m_beta * dt) * du 
                          - (m_gamma / m_beta - 1.0) * vn;
                double vp = vn + dv;
                node.set(dofV, vp);
            }
        }
    }
    
    // Update mesh
    mesh.Update(fem.GetTime());
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::UpdateIncrements(std::vector<double>& Ui, std::vector<double>& ui)
{
    for (int i = 0; i < m_neq; ++i)
    {
        Ui[i] += ui[i];
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::InternalForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
        {
            dom.InternalForces(R);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::ExternalForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
    // Boundary conditions
    for (int i = 0; i < fem.BoundaryConditions(); ++i)
    {
        FEBoundaryCondition* pbc = fem.BoundaryCondition(i);
        if (pbc && pbc->IsActive())
        {
            // pbc->LoadVector(R);
        }
    }
    
    // Model loads
    for (int i = 0; i < fem.ModelLoads(); ++i)
    {
        FEModelLoad* pml = fem.ModelLoad(i);
        if (pml && pml->IsActive())
        {
            // pml->LoadVector(R);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::ContactForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
        if (pci && pci->IsActive())
        {
            // pci->LoadVector(R);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc && plc->IsActive())
        {
            // plc->LoadVector(R, tp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::InertialForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    int NN = mesh.Nodes();
    
    // M * a where a is calculated from Newmark
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        // Get nodal mass (simplified - should come from element mass matrix)
        double mass = 1.0;  // TODO: Calculate properly
        
        for (int j = 0; j < 3; ++j)
        {
            int dofU = m_dofU[j];
            int dofV = m_dofV[j];
            int eq = node.m_ID[dofU];
            
            if (eq >= 0)
            {
                // Calculate acceleration from Newmark
                double un = node.get(dofU);
                double vn = node.get(dofV);
                double vn_old = node.get_prev(dofV);
                
                double a = (vn - vn_old) / dt;
                
                // Add inertial force (negative because it opposes acceleration)
                R[eq] -= mass * a;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::DampingForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    int NN = mesh.Nodes();
    
    // Simplified Rayleigh damping: C * v
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofV = m_dofV[j];
            int eq = node.m_ID[dofV];
            
            if (eq >= 0)
            {
                double v = node.get(dofV);
                double damp = m_dampingRatio * v;  // Simplified
                
                R[eq] -= damp;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::ContactStiffness(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
        if (pci && pci->IsActive())
        {
            // pci->StiffnessMatrix(LS);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    
    for (int i = 0; i < fem.NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = fem.NonlinearConstraint(i);
        if (plc && plc->IsActive())
        {
            // plc->StiffnessMatrix(LS, tp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::MassMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // Assemble mass matrix from domains
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (dom.IsActive())
        {
            // dom.MassMatrix(LS);
        }
    }
}

//-----------------------------------------------------------------------------
void FEImplicitDynamicSolver::DampingMatrix(FELinearSystem& LS)
{
    // Simplified implementation
    // In practice, this would assemble Rayleigh damping: C = alpha*M + beta*K
}
