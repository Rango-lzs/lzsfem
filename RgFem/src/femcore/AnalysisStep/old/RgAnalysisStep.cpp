/*********************************************************************
 * \file   RgAnalysisStep.cpp
 * \brief  Implementation of analysis step classes
 *
 * \author 
 * \date   February 2025
 *********************************************************************/

#include "RgAnalysisStep.h"
#include "femcore/FEModel.h"
#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEModelLoad.h"
#include "logger/log.h"
#include <algorithm>

//=============================================================================
// RgAnalysisStep
//=============================================================================

RgAnalysisStep::RgAnalysisStep(FEModel* fem)
    : m_fem(fem)
    , m_stepNumber(-1)
    , m_timePeriod(1.0)
    , m_initialTimeIncrement(0.1)
    , m_minTimeIncrement(1e-5)
    , m_maxTimeIncrement(0.1)
    , m_previousStep(nullptr)
    , m_activationMode(StepActivationMode::INHERITED)
    , m_solver(nullptr)
{
}

//-----------------------------------------------------------------------------
RgAnalysisStep::~RgAnalysisStep()
{
    // Note: We don't delete BCs and loads here because they might be
    // inherited by other steps. The FEModel is responsible for cleanup.
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::AddBoundaryCondition(FEBoundaryCondition* bc)
{
    if (bc)
    {
        m_boundaryConditions.push_back(bc);
    }
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::RemoveBoundaryCondition(FEBoundaryCondition* bc)
{
    auto it = std::find(m_boundaryConditions.begin(), m_boundaryConditions.end(), bc);
    if (it != m_boundaryConditions.end())
    {
        m_boundaryConditions.erase(it);
    }
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::RemoveBoundaryCondition(const std::string& name)
{
    FEBoundaryCondition* bc = FindBoundaryCondition(name);
    if (bc)
    {
        RemoveBoundaryCondition(bc);
    }
}

//-----------------------------------------------------------------------------
FEBoundaryCondition* RgAnalysisStep::FindBoundaryCondition(const std::string& name)
{
    for (auto bc : m_boundaryConditions)
    {
        if (bc->GetName() == name)
            return bc;
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::AddLoad(FEModelLoad* load)
{
    if (load)
    {
        m_loads.push_back(load);
    }
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::RemoveLoad(FEModelLoad* load)
{
    auto it = std::find(m_loads.begin(), m_loads.end(), load);
    if (it != m_loads.end())
    {
        m_loads.erase(it);
    }
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::RemoveLoad(const std::string& name)
{
    FEModelLoad* load = FindLoad(name);
    if (load)
    {
        RemoveLoad(load);
    }
}

//-----------------------------------------------------------------------------
FEModelLoad* RgAnalysisStep::FindLoad(const std::string& name)
{
    for (auto load : m_loads)
    {
        if (load->GetName() == name)
            return load;
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
void RgAnalysisStep::InheritFromPreviousStep()
{
    m_inheritedBCs.clear();
    m_inheritedLoads.clear();

    if (!m_previousStep)
        return;

    switch (m_activationMode)
    {
        case StepActivationMode::NEW:
            // Don't inherit anything
            break;

        case StepActivationMode::INHERITED:
        {
            // Inherit all active BCs from previous step
            std::vector<FEBoundaryCondition*> prevBCs = m_previousStep->GetAllActiveBCs();
            for (auto bc : prevBCs)
            {
                // Check if this BC is not redefined in current step
                if (!FindBoundaryCondition(bc->GetName()))
                {
                    m_inheritedBCs.push_back(bc);
                }
            }

            // Inherit all active loads from previous step
            std::vector<FEModelLoad*> prevLoads = m_previousStep->GetAllActiveLoads();
            for (auto load : prevLoads)
            {
                // Check if this load is not redefined in current step
                if (!FindLoad(load->GetName()))
                {
                    m_inheritedLoads.push_back(load);
                }
            }
        }
        break;

        case StepActivationMode::REPLACE:
            // Only use BCs and loads defined in this step
            break;
    }
}

//-----------------------------------------------------------------------------
std::vector<FEBoundaryCondition*> RgAnalysisStep::GetAllActiveBCs()
{
    std::vector<FEBoundaryCondition*> allBCs;
    
    // Add inherited BCs
    allBCs.insert(allBCs.end(), m_inheritedBCs.begin(), m_inheritedBCs.end());
    
    // Add BCs defined in this step
    allBCs.insert(allBCs.end(), m_boundaryConditions.begin(), m_boundaryConditions.end());
    
    return allBCs;
}

//-----------------------------------------------------------------------------
std::vector<FEModelLoad*> RgAnalysisStep::GetAllActiveLoads()
{
    std::vector<FEModelLoad*> allLoads;
    
    // Add inherited loads
    allLoads.insert(allLoads.end(), m_inheritedLoads.begin(), m_inheritedLoads.end());
    
    // Add loads defined in this step
    allLoads.insert(allLoads.end(), m_loads.begin(), m_loads.end());
    
    return allLoads;
}

//-----------------------------------------------------------------------------
bool RgAnalysisStep::Initialize()
{
    RgLog("Initializing step '%s'\n", m_name.c_str());

    // Inherit from previous step if needed
    InheritFromPreviousStep();

    // Initialize all active BCs
    std::vector<FEBoundaryCondition*> allBCs = GetAllActiveBCs();
    for (auto bc : allBCs)
    {
        if (!bc->Init())
        {
            RgLogError("Failed to initialize BC '%s' in step '%s'\n", 
                      bc->GetName().c_str(), m_name.c_str());
            return false;
        }
    }

    // Initialize all active loads
    std::vector<FEModelLoad*> allLoads = GetAllActiveLoads();
    for (auto load : allLoads)
    {
        if (!load->Init())
        {
            RgLogError("Failed to initialize load '%s' in step '%s'\n", 
                      load->GetName().c_str(), m_name.c_str());
            return false;
        }
    }

    RgLog("  Inherited BCs: %d\n", (int)m_inheritedBCs.size());
    RgLog("  New BCs: %d\n", (int)m_boundaryConditions.size());
    RgLog("  Total active BCs: %d\n", (int)allBCs.size());
    RgLog("  Inherited loads: %d\n", (int)m_inheritedLoads.size());
    RgLog("  New loads: %d\n", (int)m_loads.size());
    RgLog("  Total active loads: %d\n", (int)allLoads.size());

    return true;
}

//-----------------------------------------------------------------------------
bool RgAnalysisStep::Activate()
{
    RgLog("Activating step '%s'\n", m_name.c_str());

    // Activate all BCs
    std::vector<FEBoundaryCondition*> allBCs = GetAllActiveBCs();
    for (auto bc : allBCs)
    {
        bc->Activate();
    }

    // Activate all loads
    std::vector<FEModelLoad*> allLoads = GetAllActiveLoads();
    for (auto load : allLoads)
    {
        load->Activate();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool RgAnalysisStep::Deactivate()
{
    RgLog("Deactivating step '%s'\n", m_name.c_str());

    // Deactivate BCs defined in this step (inherited ones stay active)
    for (auto bc : m_boundaryConditions)
    {
        bc->Deactivate();
    }

    // Deactivate loads defined in this step
    for (auto load : m_loads)
    {
        load->Deactivate();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool RgAnalysisStep::Solve()
{
    RgLogError("RgAnalysisStep::Solve() must be overridden in derived class\n");
    return false;
}

//=============================================================================
// RgStaticAnalysisStep
//=============================================================================

RgStaticAnalysisStep::RgStaticAnalysisStep(FEModel* fem)
    : RgAnalysisStep(fem)
    , m_nonlinear(true)
    , m_maxIterations(25)
    , m_convergenceTol(1e-6)
{
}

//-----------------------------------------------------------------------------
RgStaticAnalysisStep::~RgStaticAnalysisStep()
{
}

//-----------------------------------------------------------------------------
bool RgStaticAnalysisStep::Initialize()
{
    RgLog("Initializing static analysis step '%s'\n", m_name.c_str());
    RgLog("  Nonlinear: %s\n", m_nonlinear ? "Yes" : "No");
    RgLog("  Max iterations: %d\n", m_maxIterations);
    RgLog("  Convergence tolerance: %g\n", m_convergenceTol);

    return RgAnalysisStep::Initialize();
}

//-----------------------------------------------------------------------------
bool RgStaticAnalysisStep::Solve()
{
    RgLog("Solving static analysis step '%s'\n", m_name.c_str());

    if (!m_solver)
    {
        RgLogError("No solver assigned to step '%s'\n", m_name.c_str());
        return false;
    }

    // Set time parameters
    double t0 = m_fem->GetCurrentTime();
    double dt = m_initialTimeIncrement;
    double tend = t0 + m_timePeriod;

    RgLog("  Time: %g -> %g\n", t0, tend);
    RgLog("  Initial dt: %g\n", dt);

    // TODO: Implement actual solving logic
    // This is a placeholder - you'll need to integrate with your solver
    
    // Example pseudo-code:
    // while (t < tend)
    // {
    //     m_fem->SetCurrentTime(t + dt);
    //     bool converged = m_solver->SolveTimeStep();
    //     if (!converged) return false;
    //     t += dt;
    // }

    RgLog("Static analysis step completed\n");
    return true;
}

//=============================================================================
// RgDynamicAnalysisStep
//=============================================================================

RgDynamicAnalysisStep::RgDynamicAnalysisStep(FEModel* fem)
    : RgAnalysisStep(fem)
    , m_explicit(false)
    , m_alpha(0.0)
    , m_beta(0.25)
    , m_gamma(0.5)
{
}

//-----------------------------------------------------------------------------
RgDynamicAnalysisStep::~RgDynamicAnalysisStep()
{
}

//-----------------------------------------------------------------------------
bool RgDynamicAnalysisStep::Initialize()
{
    RgLog("Initializing dynamic analysis step '%s'\n", m_name.c_str());
    RgLog("  Explicit: %s\n", m_explicit ? "Yes" : "No");
    if (!m_explicit)
    {
        RgLog("  Newmark parameters: alpha=%g, beta=%g, gamma=%g\n", 
              m_alpha, m_beta, m_gamma);
    }

    return RgAnalysisStep::Initialize();
}

//-----------------------------------------------------------------------------
bool RgDynamicAnalysisStep::Solve()
{
    RgLog("Solving dynamic analysis step '%s'\n", m_name.c_str());

    if (!m_solver)
    {
        RgLogError("No solver assigned to step '%s'\n", m_name.c_str());
        return false;
    }

    // Set time parameters
    double t0 = m_fem->GetCurrentTime();
    double dt = m_initialTimeIncrement;
    double tend = t0 + m_timePeriod;

    RgLog("  Time: %g -> %g\n", t0, tend);
    RgLog("  Initial dt: %g\n", dt);

    // TODO: Implement actual solving logic for dynamic analysis
    
    RgLog("Dynamic analysis step completed\n");
    return true;
}
