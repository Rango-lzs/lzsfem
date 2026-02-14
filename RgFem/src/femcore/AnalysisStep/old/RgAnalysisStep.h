/*********************************************************************
 * \file   RgAnalysisStep.h
 * \brief  Analysis step class for managing BCs and loads per step
 *
 * \author 
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"
#include <string>
#include <vector>
#include <memory>

class FEModel;
class FEBoundaryCondition;
class FEModelLoad;
class FESolver;

/**
 * @brief Step activation mode for BCs and loads
 */
enum class StepActivationMode
{
    NEW,         // Only new BCs/loads defined in this step
    INHERITED,   // Inherit from previous step + new ones
    REPLACE      // Replace all with new definitions
};

/**
 * @brief Analysis step base class
 * 
 * This class manages boundary conditions and loads for a specific analysis step.
 * BCs and loads can be:
 * - Newly defined in this step
 * - Inherited from previous steps
 * - Modified or removed in this step
 */
class FEM_EXPORT RgAnalysisStep
{
public:
    RgAnalysisStep(FEModel* fem);
    virtual ~RgAnalysisStep();

    // Basic step properties
    void SetName(const std::string& name) { m_name = name; }
    std::string GetName() const { return m_name; }
    
    void SetStepNumber(int n) { m_stepNumber = n; }
    int GetStepNumber() const { return m_stepNumber; }

    // Time control
    void SetTimePeriod(double t) { m_timePeriod = t; }
    double GetTimePeriod() const { return m_timePeriod; }
    
    void SetInitialTimeIncrement(double dt) { m_initialTimeIncrement = dt; }
    double GetInitialTimeIncrement() const { return m_initialTimeIncrement; }
    
    void SetMinTimeIncrement(double dt) { m_minTimeIncrement = dt; }
    double GetMinTimeIncrement() const { return m_minTimeIncrement; }
    
    void SetMaxTimeIncrement(double dt) { m_maxTimeIncrement = dt; }
    double GetMaxTimeIncrement() const { return m_maxTimeIncrement; }

    // Boundary condition management
    void AddBoundaryCondition(FEBoundaryCondition* bc);
    void RemoveBoundaryCondition(FEBoundaryCondition* bc);
    void RemoveBoundaryCondition(const std::string& name);
    int BoundaryConditions() const { return (int)m_boundaryConditions.size(); }
    FEBoundaryCondition* GetBoundaryCondition(int i) { return m_boundaryConditions[i]; }
    const FEBoundaryCondition* GetBoundaryCondition(int i) const { return m_boundaryConditions[i]; }
    FEBoundaryCondition* FindBoundaryCondition(const std::string& name);

    // Load management
    void AddLoad(FEModelLoad* load);
    void RemoveLoad(FEModelLoad* load);
    void RemoveLoad(const std::string& name);
    int Loads() const { return (int)m_loads.size(); }
    FEModelLoad* GetLoad(int i) { return m_loads[i]; }
    const FEModelLoad* GetLoad(int i) const { return m_loads[i]; }
    FEModelLoad* FindLoad(const std::string& name);

    // Inheritance management
    void SetActivationMode(StepActivationMode mode) { m_activationMode = mode; }
    StepActivationMode GetActivationMode() const { return m_activationMode; }
    
    void SetPreviousStep(RgAnalysisStep* prevStep) { m_previousStep = prevStep; }
    RgAnalysisStep* GetPreviousStep() { return m_previousStep; }
    
    // Inherit BCs and loads from previous step
    void InheritFromPreviousStep();
    
    // Get all active BCs (including inherited)
    std::vector<FEBoundaryCondition*> GetAllActiveBCs();
    
    // Get all active loads (including inherited)
    std::vector<FEModelLoad*> GetAllActiveLoads();

    // Step lifecycle
    virtual bool Initialize();
    virtual bool Activate();
    virtual bool Deactivate();
    virtual bool Solve();

    // Model access
    FEModel* GetFEModel() { return m_fem; }
    const FEModel* GetFEModel() const { return m_fem; }

    // Solver management
    void SetSolver(FESolver* solver) { m_solver = solver; }
    FESolver* GetSolver() { return m_solver; }

protected:
    FEModel* m_fem;
    std::string m_name;
    int m_stepNumber;

    // Time control parameters
    double m_timePeriod;
    double m_initialTimeIncrement;
    double m_minTimeIncrement;
    double m_maxTimeIncrement;

    // BCs and loads defined in this step
    std::vector<FEBoundaryCondition*> m_boundaryConditions;
    std::vector<FEModelLoad*> m_loads;

    // Inherited BCs and loads (references, not owned)
    std::vector<FEBoundaryCondition*> m_inheritedBCs;
    std::vector<FEModelLoad*> m_inheritedLoads;

    // Step relationship
    RgAnalysisStep* m_previousStep;
    StepActivationMode m_activationMode;

    // Solver
    FESolver* m_solver;
};

/**
 * @brief Static analysis step
 */
class FEM_EXPORT RgStaticAnalysisStep : public RgAnalysisStep
{
public:
    RgStaticAnalysisStep(FEModel* fem);
    virtual ~RgStaticAnalysisStep();

    bool Initialize() override;
    bool Solve() override;

    // Static analysis specific parameters
    void SetNonlinear(bool b) { m_nonlinear = b; }
    bool IsNonlinear() const { return m_nonlinear; }

    void SetMaxIterations(int n) { m_maxIterations = n; }
    int GetMaxIterations() const { return m_maxIterations; }

    void SetConvergenceTolerance(double tol) { m_convergenceTol = tol; }
    double GetConvergenceTolerance() const { return m_convergenceTol; }

private:
    bool m_nonlinear;
    int m_maxIterations;
    double m_convergenceTol;
};

/**
 * @brief Dynamic analysis step
 */
class FEM_EXPORT RgDynamicAnalysisStep : public RgAnalysisStep
{
public:
    RgDynamicAnalysisStep(FEModel* fem);
    virtual ~RgDynamicAnalysisStep();

    bool Initialize() override;
    bool Solve() override;

    // Dynamic analysis specific parameters
    void SetExplicit(bool b) { m_explicit = b; }
    bool IsExplicit() const { return m_explicit; }

    void SetAlpha(double alpha) { m_alpha = alpha; }
    double GetAlpha() const { return m_alpha; }

    void SetBeta(double beta) { m_beta = beta; }
    double GetBeta() const { return m_beta; }

    void SetGamma(double gamma) { m_gamma = gamma; }
    double GetGamma() const { return m_gamma; }

private:
    bool m_explicit;
    double m_alpha;  // Newmark parameters
    double m_beta;
    double m_gamma;
};
