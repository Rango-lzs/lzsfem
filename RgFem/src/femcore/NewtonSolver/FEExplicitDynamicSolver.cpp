/*********************************************************************
 * \file   FEExplicitDynamicSolver.cpp
 * \brief  Implementation of explicit dynamic solver
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "FEExplicitDynamicSolver.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/Domain/RgDomain.h"
#include "femcore/FEAnalysis/FEAnalysis.h"
#include "logger/log.h"
#include <cmath>

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEExplicitDynamicSolver, FESolver)
    ADD_PARAMETER(m_lumpedMass, "lumped_mass");
    ADD_PARAMETER(m_autoTimeStep, "auto_time_step");
    ADD_PARAMETER(m_timeStepScale, "time_step_scale");
    ADD_PARAMETER(m_maxTimeStep, "max_time_step");
    ADD_PARAMETER(m_checkStability, "check_stability");
    ADD_PARAMETER(m_useDamping, "use_damping");
    ADD_PARAMETER(m_dampingCoeff, "damping_coeff");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEExplicitDynamicSolver::FEExplicitDynamicSolver()
{
    m_scheme = CENTRAL_DIFFERENCE;
    m_lumpedMass = true;
    m_autoTimeStep = false;
    m_timeStepScale = 0.9;  // 90% of critical for safety
    m_maxTimeStep = 1e10;
    m_checkStability = true;
    m_useDamping = false;
    m_dampingCoeff = 0.0;
    
    m_pM = nullptr;
    m_dtCritical = 0.0;
    m_initialized = false;
}

//-----------------------------------------------------------------------------
FEExplicitDynamicSolver::~FEExplicitDynamicSolver()
{
    if (m_pM) delete m_pM;
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::Serialize(DumpStream& ar)
{
    FESolver::Serialize(ar);
    
    ar & m_Un & m_Vn & m_An & m_Vnhalf;
    ar & m_Fint & m_Fext & m_Fcon;
    ar & m_Mi;
    ar & m_dtCritical;
    ar & m_initialized;
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::Init()
{
    if (!FESolver::Init()) return false;
    
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
    
    m_dofA.Clear();
    m_dofA.AddDof(dofs.GetDOF("ax"));
    m_dofA.AddDof(dofs.GetDOF("ay"));
    m_dofA.AddDof(dofs.GetDOF("az"));
    
    // Allocate state vectors
    int neq = m_neq;
    m_Un.assign(neq, 0.0);
    m_Vn.assign(neq, 0.0);
    m_An.assign(neq, 0.0);
    m_Vnhalf.assign(neq, 0.0);
    
    m_Fint.assign(neq, 0.0);
    m_Fext.assign(neq, 0.0);
    m_Fcon.assign(neq, 0.0);
    
    // Calculate mass matrix
    if (!CalculateMassMatrix())
        return false;
    
    // Calculate critical time step
    m_dtCritical = CalculateCriticalTimeStep();
    RgLog("Critical time step: %.6e\n", m_dtCritical);
    
    m_initialized = true;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::InitEquations()
{
    return FESolver::InitEquations();
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::InitStep(double time)
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    
    tp.currentTime = time;
    
    if (!m_initialized)
    {
        InitializeVelocity();
        m_initialized = true;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::InitializeVelocity()
{
    // For central difference, we need v(n-1/2)
    // We can initialize it using: v(n-1/2) = v(n) - 0.5*dt*a(n)
    
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    int NN = mesh.Nodes();
    
    // First calculate initial accelerations
    CalculateAccelerations();
    
    // Then calculate v(n-1/2)
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofV = m_dofV[j];
            int dofA = m_dofA[j];
            int eq = node.m_ID[dofV];
            
            if (eq >= 0)
            {
                double vn = node.get(dofV);
                double an = node.get(dofA);
                
                m_Vnhalf[eq] = vn - 0.5 * dt * an;
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::SolveStep()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Check stability if enabled
    if (m_checkStability && !CheckStability())
    {
        RgLogError("Time step exceeds stability limit");
        return false;
    }
    
    // Adjust time step if auto mode
    if (m_autoTimeStep)
    {
        double dt_safe = m_timeStepScale * m_dtCritical;
        if (dt_safe < m_maxTimeStep)
        {
            step->m_dt = dt_safe;
            dt = dt_safe;
        }
    }
    
    // Perform time integration based on scheme
    switch (m_scheme)
    {
        case CENTRAL_DIFFERENCE:
            CentralDifferenceStep();
            break;
            
        case FORWARD_EULER:
            ForwardEulerStep();
            break;
            
        case VELOCITY_VERLET:
            VelocityVerletStep();
            break;
    }
    
    // Update model
    fem.Update();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::CentralDifferenceStep()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Step 1: Calculate accelerations at time n
    CalculateAccelerations();
    
    // Step 2: Update velocities to n+1/2
    UpdateVelocities();
    
    // Step 3: Update displacements to n+1
    UpdateDisplacements();
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::ForwardEulerStep()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Calculate accelerations
    CalculateAccelerations();
    
    // Update velocities: v(n+1) = v(n) + dt*a(n)
    for (int i = 0; i < m_neq; ++i)
    {
        m_Vn[i] += dt * m_An[i];
    }
    
    // Update displacements: u(n+1) = u(n) + dt*v(n+1)
    for (int i = 0; i < m_neq; ++i)
    {
        m_Un[i] += dt * m_Vn[i];
    }
    
    // Update model
    std::vector<double> du(m_neq);
    for (int i = 0; i < m_neq; ++i)
    {
        du[i] = dt * m_Vn[i];
    }
    Update(du);
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::VelocityVerletStep()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Step 1: Update positions
    // u(n+1) = u(n) + v(n)*dt + 0.5*a(n)*dt^2
    for (int i = 0; i < m_neq; ++i)
    {
        m_Un[i] += m_Vn[i] * dt + 0.5 * m_An[i] * dt * dt;
    }
    
    std::vector<double> du(m_neq);
    for (int i = 0; i < m_neq; ++i)
    {
        du[i] = m_Vn[i] * dt + 0.5 * m_An[i] * dt * dt;
    }
    Update(du);
    
    // Step 2: Calculate new accelerations
    std::vector<double> an_old = m_An;
    CalculateAccelerations();
    
    // Step 3: Update velocities
    // v(n+1) = v(n) + 0.5*(a(n) + a(n+1))*dt
    for (int i = 0; i < m_neq; ++i)
    {
        m_Vn[i] += 0.5 * (an_old[i] + m_An[i]) * dt;
    }
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::CalculateAccelerations()
{
    FEModel& fem = *GetFEModel();
    
    // Zero force vectors
    std::fill(m_Fint.begin(), m_Fint.end(), 0.0);
    std::fill(m_Fext.begin(), m_Fext.end(), 0.0);
    std::fill(m_Fcon.begin(), m_Fcon.end(), 0.0);
    
    // Create global vectors
    FEGlobalVector RHS_int(fem, m_Fint, m_Fd);
    FEGlobalVector RHS_ext(fem, m_Fext, m_Fd);
    FEGlobalVector RHS_con(fem, m_Fcon, m_Fd);
    
    // Calculate forces
    InternalForces(RHS_int);
    ExternalForces(RHS_ext);
    ContactForces(RHS_con);
    
    // Calculate accelerations: a = M^-1 * (Fext - Fint + Fcon)
    for (int i = 0; i < m_neq; ++i)
    {
        double F_total = m_Fext[i] - m_Fint[i] + m_Fcon[i];
        
        // Apply damping if enabled
        if (m_useDamping)
        {
            F_total -= m_dampingCoeff * m_Vn[i];
        }
        
        m_An[i] = F_total / m_Mi[i];
    }
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::UpdateVelocities()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    // Update velocities: v(n+1/2) = v(n-1/2) + a(n)*dt
    for (int i = 0; i < m_neq; ++i)
    {
        m_Vnhalf[i] += m_An[i] * dt;
    }
    
    // Also update nodal velocities (for output)
    FEMesh& mesh = fem.GetMesh();
    int NN = mesh.Nodes();
    
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofV = m_dofV[j];
            int eq = node.m_ID[dofV];
            
            if (eq >= 0)
            {
                node.set(dofV, m_Vnhalf[eq]);
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::UpdateDisplacements()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    int NN = mesh.Nodes();
    
    // Update displacements: u(n+1) = u(n) + v(n+1/2)*dt
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        
        for (int j = 0; j < 3; ++j)
        {
            int dofU = m_dofU[j];
            int eq = node.m_ID[dofU];
            
            if (eq >= 0)
            {
                double du = m_Vnhalf[eq] * dt;
                m_Un[eq] += du;
                
                node.set(dofU, m_Un[eq]);
                
                // Update position
                if (j == 0) node.m_rt.x = node.m_r0.x + m_Un[eq];
                if (j == 1) node.m_rt.y = node.m_r0.y + m_Un[eq];
                if (j == 2) node.m_rt.z = node.m_r0.z + m_Un[eq];
            }
        }
    }
    
    // Update mesh
    mesh.Update(fem.GetTime());
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::CheckStability()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* step = fem.GetCurrentStep();
    
    double dt = step->m_dt;
    
    if (dt > m_dtCritical)
    {
        RgLogError("Time step %.6e exceeds critical value %.6e", dt, m_dtCritical);
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
double FEExplicitDynamicSolver::CalculateCriticalTimeStep()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    double dt_min = 1e20;
    
    // Calculate critical time step for each element
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        if (!dom.IsActive()) continue;
        
        for (int j = 0; j < dom.Elements(); ++j)
        {
            RgElement& el = dom.ElementRef(j);
            
            // Get element characteristic length
            double h = 1.0;  // TODO: Calculate actual element size
            
            // Get material wave speed
            double c = 1.0;  // TODO: Get from material (sqrt(E/rho))
            
            // CFL condition: dt <= h / c
            double dt_elem = h / c;
            
            if (dt_elem < dt_min)
                dt_min = dt_elem;
        }
    }
    
    return dt_min;
}

//-----------------------------------------------------------------------------
bool FEExplicitDynamicSolver::CalculateMassMatrix()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    if (m_lumpedMass)
    {
        // Lumped mass matrix (diagonal)
        m_Mi.assign(m_neq, 0.0);
        
        // Assemble from domains
        for (int i = 0; i < mesh.Domains(); ++i)
        {
            RgDomain& dom = mesh.Domain(i);
            if (dom.IsActive())
            {
                // dom.LumpedMassMatrix(m_Mi);
                
                // Simplified: uniform mass
                for (int j = 0; j < dom.Elements(); ++j)
                {
                    RgElement& el = dom.ElementRef(j);
                    int nen = el.NodeSize();
                    double mass = 1.0 / nen;  // Distribute uniformly
                    
                    for (int k = 0; k < nen; ++k)
                    {
                        FENode& node = mesh.Node(el.getNodeId(k));
                        
                        for (int l = 0; l < 3; ++l)
                        {
                            int eq = node.m_ID[m_dofU[l]];
                            if (eq >= 0)
                                m_Mi[eq] += mass;
                        }
                    }
                }
            }
        }
        
        // Check for zero mass
        for (int i = 0; i < m_neq; ++i)
        {
            if (m_Mi[i] < 1e-20)
            {
                RgLogError("Zero mass detected at equation %d", i);
                return false;
            }
        }
    }
    else
    {
        // Consistent mass matrix (full)
        // Create matrix
        FECoreKernel& fecore = FECoreKernel::GetInstance();
        m_pM = fecore.CreateGlobalMatrix(MatType());
        
        if (!m_pM)
        {
            RgLogError("Failed to create mass matrix");
            return false;
        }
        
        // Build profile
        BuildMatrixProfile(*m_pM, true);
        
        // Assemble mass matrix
        // TODO: Implement consistent mass assembly
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEExplicitDynamicSolver::InternalForces(FEGlobalVector& R)
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
void FEExplicitDynamicSolver::ExternalForces(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    
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
void FEExplicitDynamicSolver::ContactForces(FEGlobalVector& R)
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
void FEExplicitDynamicSolver::Update(std::vector<double>& u)
{
    // Already updated in UpdateDisplacements
}

//-----------------------------------------------------------------------------
FEGlobalMatrix* FEExplicitDynamicSolver::GetStiffnessMatrix()
{
    // Explicit solver doesn't use stiffness matrix
    return nullptr;
}

//-----------------------------------------------------------------------------
std::vector<double> FEExplicitDynamicSolver::GetLoadVector()
{
    return m_Fext;
}

//-----------------------------------------------------------------------------
LinearSolver* FEExplicitDynamicSolver::GetLinearSolver()
{
    // Explicit solver doesn't use linear solver
    return nullptr;
}
