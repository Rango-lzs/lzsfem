/*********************************************************************
 * \file   SolverUsageExample.cpp
 * \brief  Examples of how to use the refactored solver system
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "FEStaticSolver.h"
#include "FEImplicitDynamicSolver.h"
#include "FEExplicitDynamicSolver.h"
#include "femcore/FEModel.h"
#include "femcore/FEAnalysis/FEAnalysis.h"

//=============================================================================
// Example 1: Static Analysis
//=============================================================================
void Example_StaticAnalysis()
{
    // Create FE model
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (load mesh, materials, boundary conditions, etc.)
    
    // Create analysis step
    FEAnalysis* step = new FEAnalysis(&model);
    step->SetName("Static Analysis");
    
    // Create static solver
    FEStaticSolver* solver = new FEStaticSolver();
    
    // Configure solver parameters
    solver->m_Rtol = 1e-5;        // Residual tolerance
    solver->m_Etol = 1e-4;        // Energy tolerance
    solver->m_Dtol = 1e-3;        // Displacement tolerance
    solver->m_maxref = 15;        // Max reformations per step
    solver->m_breformtimestep = true;  // Reform at start of step
    
    // Set linear solver
    // solver->m_plinsolve = new PardisoSolver();
    
    // Attach solver to step
    step->SetFESolver(solver);
    
    // Add step to model
    model.AddStep(step);
    
    // Initialize and solve
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 2: Implicit Dynamic Analysis (Newmark)
//=============================================================================
void Example_ImplicitDynamic_Newmark()
{
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (setup model)
    
    FEAnalysis* step = new FEAnalysis(&model);
    step->SetName("Impact Analysis");
    
    // Create implicit dynamic solver
    FEImplicitDynamicSolver* solver = new FEImplicitDynamicSolver();
    
    // Configure time integration
    solver->m_timeScheme = NEWMARK;
    solver->m_beta = 0.25;        // Average acceleration
    solver->m_gamma = 0.5;
    
    // Configure convergence
    solver->m_Rtol = 1e-5;
    solver->m_Etol = 1e-4;
    solver->m_Dtol = 1e-3;
    solver->m_maxref = 15;
    
    // Enable damping (optional)
    solver->m_useDamping = true;
    solver->m_dampingRatio = 0.05;  // 5% damping
    
    // Attach to step
    step->SetFESolver(solver);
    model.AddStep(step);
    
    // Initialize and solve
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 3: Implicit Dynamic Analysis (Generalized-Alpha)
//=============================================================================
void Example_ImplicitDynamic_GenAlpha()
{
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (setup model)
    
    FEAnalysis* step = new FEAnalysis(&model);
    step->SetName("Vibration Analysis");
    
    // Create implicit dynamic solver
    FEImplicitDynamicSolver* solver = new FEImplicitDynamicSolver();
    
    // Configure time integration
    solver->m_timeScheme = GENERALIZED_ALPHA;
    solver->m_rhoi = 0.8;  // Spectral radius (controls high-frequency damping)
    // Note: alphaf and alpham are calculated automatically from rhoi
    
    // Configure convergence
    solver->m_Rtol = 1e-5;
    solver->m_Dtol = 1e-3;
    solver->m_maxref = 15;
    
    // Attach to step
    step->SetFESolver(solver);
    model.AddStep(step);
    
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 4: Explicit Dynamic Analysis
//=============================================================================
void Example_ExplicitDynamic()
{
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (setup model)
    
    FEAnalysis* step = new FEAnalysis(&model);
    step->SetName("Wave Propagation");
    
    // Create explicit dynamic solver
    FEExplicitDynamicSolver* solver = new FEExplicitDynamicSolver();
    
    // Configure time integration
    solver->m_scheme = CENTRAL_DIFFERENCE;
    solver->m_lumpedMass = true;   // Use lumped mass (faster)
    
    // Configure automatic time stepping
    solver->m_autoTimeStep = true;
    solver->m_timeStepScale = 0.9; // Safety factor (90% of critical)
    solver->m_checkStability = true;
    
    // Attach to step
    step->SetFESolver(solver);
    model.AddStep(step);
    
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 5: Multi-Step Analysis
//=============================================================================
void Example_MultiStep()
{
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (setup model)
    
    // Step 1: Static preload
    {
        FEAnalysis* step1 = new FEAnalysis(&model);
        step1->SetName("Preload");
        step1->m_ntime = 10;      // 10 time steps
        step1->m_dt = 0.1;        // Time increment
        
        FEStaticSolver* solver = new FEStaticSolver();
        solver->m_Rtol = 1e-5;
        solver->m_Dtol = 1e-3;
        
        step1->SetFESolver(solver);
        model.AddStep(step1);
    }
    
    // Step 2: Dynamic impact
    {
        FEAnalysis* step2 = new FEAnalysis(&model);
        step2->SetName("Impact");
        step2->m_ntime = 1000;    // 1000 time steps
        step2->m_dt = 1e-6;       // Small time step for impact
        
        FEExplicitDynamicSolver* solver = new FEExplicitDynamicSolver();
        solver->m_scheme = CENTRAL_DIFFERENCE;
        solver->m_lumpedMass = true;
        solver->m_autoTimeStep = true;
        
        step2->SetFESolver(solver);
        model.AddStep(step2);
    }
    
    // Step 3: Relaxation
    {
        FEAnalysis* step3 = new FEAnalysis(&model);
        step3->SetName("Relaxation");
        step3->m_ntime = 100;
        step3->m_dt = 0.01;
        
        FEImplicitDynamicSolver* solver = new FEImplicitDynamicSolver();
        solver->m_timeScheme = GENERALIZED_ALPHA;
        solver->m_rhoi = 0.8;
        solver->m_useDamping = true;
        solver->m_dampingRatio = 0.1;
        
        step3->SetFESolver(solver);
        model.AddStep(step3);
    }
    
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 6: Solver with Line Search and QN Strategy
//=============================================================================
void Example_AdvancedSolverSettings()
{
    FEModel model;
    model.SetActiveModule("solid");
    
    // ... (setup model)
    
    FEAnalysis* step = new FEAnalysis(&model);
    
    FEStaticSolver* solver = new FEStaticSolver();
    
    // Configure Newton iteration
    solver->m_maxref = 15;
    solver->m_Rtol = 1e-5;
    solver->m_Etol = 1e-4;
    solver->m_Dtol = 1e-3;
    
    // Configure reformation strategy
    solver->m_breformtimestep = true;   // Reform at start of step
    solver->m_bdivreform = true;        // Reform when diverging
    solver->m_breformAugment = false;   // Don't reform after augmentation
    
    // Set QN strategy
    solver->SetDefaultStrategy(QN_BFGS);  // Use BFGS updates
    
    // Configure line search
    FELineSearch* ls = solver->GetLineSearch();
    if (ls)
    {
        ls->m_maxiter = 10;      // Max line search iterations
        ls->m_tolerance = 0.9;   // Line search tolerance
    }
    
    // Check for zero diagonals
    solver->CheckZeroDiagonal(true, 1e-10);
    
    step->SetFESolver(solver);
    model.AddStep(step);
    
    if (model.Init())
    {
        model.Solve();
    }
}

//=============================================================================
// Example 7: Accessing Solution Data
//=============================================================================
void Example_AccessSolutionData()
{
    FEModel model;
    // ... (setup and solve)
    
    // Get current step
    FEAnalysis* step = model.GetCurrentStep();
    FESolver* solver = step->GetFESolver();
    
    // For static or implicit dynamic solver
    if (FENewtonSolver* newton = dynamic_cast<FENewtonSolver*>(solver))
    {
        // Get solution vector
        std::vector<double> U = newton->GetSolutionVector();
        
        // Get load vector
        std::vector<double> F = newton->GetLoadVector();
        
        // Get stiffness matrix
        FEGlobalMatrix* K = newton->GetStiffnessMatrix();
        
        // Check if converged
        // ... (convergence info available in solver)
    }
    
    // For explicit dynamic solver
    if (FEExplicitDynamicSolver* expl = dynamic_cast<FEExplicitDynamicSolver*>(solver))
    {
        // Get critical time step
        // double dt_crit = expl->m_dtCritical;
    }
    
    // Access nodal displacements
    FEMesh& mesh = model.GetMesh();
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        Vector3d u = node.m_rt - node.m_r0;  // Displacement
        // ... use displacement
    }
}
