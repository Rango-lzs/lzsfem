#pragma once
#include "femcore/FEObjectBase.h"

#include <string>
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class RgDomain;
class DumpStream;
class FEStepComponent;
class FETimeStepController;
class RgBoundaryCondition;
class RgLoad;

//-----------------------------------------------------------------------------
//! Step activation mode for BCs and loads
enum class StepActivationMode
{
    NEW,        // Only new BCs/loads defined in this step
    INHERITED,  // Inherit from previous step + new ones (default)
    REPLACE     // Replace all with new definitions
};

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
//! Extended to support BC and Load management per step
class FEM_EXPORT FEAnalysis : public FEObjectBase
{
    DECLARE_META_CLASS(FEAnalysis, FEObjectBase);

public:
    //! constructor
    FEAnalysis();

    //! destructor
    virtual ~FEAnalysis() = 0;

    //! Initialization
    virtual bool Init() override;

    //! Reset analysis data
    virtual void Reset()
    {
    }

    //! Serialize data from and to a binary archive
    virtual void Serialize(DumpStream& ar) override;

    //! copy data from another analysis
    void CopyFrom(FEAnalysis* step);

    //-----------------------------------------------------------------------------
    //! See if this step is active
    bool IsActive();

    //-----------------------------------------------------------------------------
    //! This function gets called right before the step needs to be solved.
    bool Activate();

public:  // Solver management
    void SetFESolver(FESolver* psolver);
    FESolver* GetFESolver();

    // initialize the solver
    bool InitSolver();

    virtual bool Solve();

    // Call the FE Solver to solve the time step
    // Returns an error code
    // 0 = all is well, continue
    // 1 = solver has failed, but try auto-time step
    // 2 = abort
    int SolveTimeStep();

public:  // Domain management
    //! Get active domains
    int Domains()
    {
        return (int)m_Dom.size();
    }

    //! Get active domain
    RgDomain* Domain(int i);

    //! Add a domain
    void AddDomain(int i)
    {
        m_Dom.push_back(i);
    }

    //! clear all domains
    void ClearDomains()
    {
        m_Dom.clear();
    }

public:  // Step component management
    //! add a step component
    void AddStepComponent(FEStepComponent* pmc);

    //! return number of model components
    int StepComponents() const;

    //! get a step component
    FEStepComponent* GetStepComponent(int i);

public:  // NEW: Boundary condition management
    //! Add a boundary condition to this step
    void AddBoundaryCondition(RgBoundaryCondition* bc);

    //! Remove a boundary condition from this step
    void RemoveBoundaryCondition(RgBoundaryCondition* bc);
    void RemoveBoundaryCondition(const std::string& name);

    //! Get number of BCs defined in this step (not including inherited)
    int BoundaryConditions() const
    {
        return (int)m_BC.size();
    }

    //! Get a BC defined in this step
    RgBoundaryCondition* GetBoundaryCondition(int i)
    {
        return m_BC[i];
    }
    const RgBoundaryCondition* GetBoundaryCondition(int i) const
    {
        return m_BC[i];
    }

    //! Find a BC by name
    RgBoundaryCondition* FindBoundaryCondition(const std::string& name);

    //! Get all active BCs (including inherited)
    std::vector<RgBoundaryCondition*> GetAllActiveBCs();

public:  // NEW: Load management
    //! Add a load to this step
    void AddLoad(RgLoad* load);

    //! Remove a load from this step
    void RemoveLoad(RgLoad* load);
    void RemoveLoad(const std::string& name);

    //! Get number of loads defined in this step (not including inherited)
    int Loads() const
    {
        return (int)m_Loads.size();
    }

    //! Get a load defined in this step
    RgLoad* GetLoad(int i)
    {
        return m_Loads[i];
    }
    const RgLoad* GetLoad(int i) const
    {
        return m_Loads[i];
    }

    //! Find a load by name
    RgLoad* FindLoad(const std::string& name);

    //! Get all active loads (including inherited)
    std::vector<RgLoad*> GetAllActiveLoads();

public:  // NEW: Inheritance management
    //! Set the activation mode for this step
    void SetActivationMode(StepActivationMode mode)
    {
        m_activationMode = mode;
    }
    StepActivationMode GetActivationMode() const
    {
        return m_activationMode;
    }

    //! Set the previous step (for inheritance)
    void SetPreviousStep(FEAnalysis* prevStep)
    {
        m_previousStep = prevStep;
    }
    FEAnalysis* GetPreviousStep()
    {
        return m_previousStep;
    }
    const FEAnalysis* GetPreviousStep() const
    {
        return m_previousStep;
    }

    //! Inherit BCs and loads from previous step
    void InheritFromPreviousStep();

    //! Activate all BCs and loads for this step
    bool ActivateBCsAndLoads();

    //! Deactivate BCs and loads defined in this step
    void DeactivateBCsAndLoads();

public:  // Plot and output settings
    //! sets the plot level
    void SetPlotLevel(int n);

    //! sets the plot stride
    void SetPlotStride(int n);

    //! sets the plot range
    void SetPlotRange(int n0, int n1);

    //! sets the zero-state plot flag
    void SetPlotZeroState(bool b);

    //! sets the plot hint
    void SetPlotHint(int plotHint);

    //! get the plot hint
    int GetPlotHint() const;

    //! get the plot level
    int GetPlotLevel();

    //! Set the output level
    void SetOutputLevel(int n);

    //! Get the output level
    int GetOutputLevel();

public:
    // --- Control Data ---
    //{
    int m_nanalysis;         //!< analysis type
    bool m_badaptorReSolve;  //!< resolve analysis after mesh adaptor phase
                             //}

                             // --- Time Step Data ---
    //{
    int m_ntime;          //!< number of time steps
    double m_final_time;  //!< end time for this time step
    double m_dt;          //!< current time step
    double m_dt0;         //!< initial time step size
    double m_tstart;      //!< start time
    double m_tend;        //!< end time

    FETimeStepController* m_timeController = nullptr;
    //}

    // --- Quasi-Newton Solver Variables ---
    //{
    int m_ntotrhs;     //!< total nr of right hand side evaluations
    int m_ntotref;     //!< total nr of stiffness reformations
    int m_ntotiter;    //!< total nr of non-linear iterations
    int m_ntimesteps;  //!< number of completed time steps
                       //}

                       // --- I/O Data ---
    //{
    int m_nplot;          //!< plot level
    int m_noutput;        //!< data output level
    int m_nplot_stride;   //!< stride for plotting
    int m_nplotRange[2];  //!< plot range
    bool m_bplotZero;     //!< Force plotting of time step "zero"
    int m_plotHint;       //!< the plot mode
    //}

private:
    // the FE solver
    FESolver* m_psolver = nullptr;  //!< pointer to solver class that will solve this step.
    bool m_bactive;                 //!< activation flag

protected:
    std::vector<int> m_Dom;              //!< list of active domains for this analysis
    std::vector<FEStepComponent*> m_MC;  //!< array of model components active during this step

    // NEW: BC and Load management
    std::vector<RgBoundaryCondition*> m_BC;           //!< BCs defined in this step
    std::vector<RgLoad*> m_Loads;                //!< Loads defined in this step
    std::vector<RgBoundaryCondition*> m_inheritedBC;  //!< BCs inherited from previous step
    std::vector<RgLoad*> m_inheritedLoads;       //!< Loads inherited from previous step

    // NEW: Step relationship
    FEAnalysis* m_previousStep;           //!< pointer to previous step
    StepActivationMode m_activationMode;  //!< how to handle BC/Load inheritance

    DECLARE_PARAM_LIST();
};