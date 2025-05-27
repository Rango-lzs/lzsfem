#pragma once
#include "femcore/FEObjectBase.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FESolver;
class FEDomain;
class DumpStream;
class FEStepComponent;
class FETimeStepController;

//-----------------------------------------------------------------------------
//! Base class for finite element analysis
class FEM_EXPORT FEAnalysis : public FEObjectBase
{
	DECLARE_META_CLASS(FEAnalysis, FEObjectBase);

public:
	//! constructor
	FEAnalysis(FEModel* pfem);

	//! destructor
	virtual ~FEAnalysis() =0;

	//! Initialization
	virtual bool Init() override;

	//! Reset analysis data
    virtual void Reset(){}

	//! Serialize data from and to a binary archive
	virtual void Serialize(DumpStream& ar) override;

	//! copy data from another analysis
	void CopyFrom(FEAnalysis* step);

public:
	void SetFESolver(FESolver* psolver);

	FESolver* GetFESolver();

	// initialize the solver
    bool InitSolver();

	bool virtual Solve();

    // Call the FE Solver to solve the time step
    // Returns an error code
    // 0 = all is well, continue
    // 1 = solver has failed, but try auto-time step
    // 2 = abort
    int SolveTimeStep();

public:
	//! Get active domains
	int Domains() { return (int)m_Dom.size(); }

	//! Get active domain
	FEDomain* Domain(int i);

	//! Add a domain
	void AddDomain(int i) { m_Dom.push_back(i); }

	//! clear all domains
	void ClearDomains() { m_Dom.clear(); }

public:
	//! add a step component
	void AddStepComponent(FEStepComponent* pmc);

	//! return number of model components
	int StepComponents() const;

	//! get a step component
	FEStepComponent* GetStepComponent(int i);

public:
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
		int		m_nanalysis;		//!< analysis type
		bool	m_badaptorReSolve;	//!< resolve analysis after mesh adaptor phase
	//}

	// --- Time Step Data ---
	//{
		int		m_ntime;		//时间步数
		double	m_final_time;	//!< end time for this time step
		double	m_dt;			//!< current time step 
		double	m_dt0;			//!< initial time step size
		double	m_tstart;		//!< start time
		double	m_tend;			//!< end time

		FETimeStepController* m_timeController;
	//}

	// --- Quasi-Newton Solver Variables ---
	//{
		int		m_ntotrhs;		//!< total nr of right hand side evaluations
		int		m_ntotref;		//!< total nr of stiffness reformations
		int		m_ntotiter;		//!< total nr of non-linear iterations
		int		m_ntimesteps;	//时间增量步数
	//}

	// --- I/O Data ---
	//{
		int		m_nplot;		//!< plot level
		int		m_noutput;		//!< data output level
		int		m_nplot_stride;	//!< stride for plotting
		int		m_nplotRange[2];	//!< plot range
		bool	m_bplotZero;		//!< Force plotting of time step "zero"
		int		m_plotHint;			//!< the plot mode
	//}

private:
	// the FE solver
	FESolver*	m_psolver;	//!< pointer to solver class that will solve this step.
	bool		m_bactive;	//!< activation flag

protected:
	std::vector<int>				m_Dom;	//!< list of active domains for this analysis
	std::vector<FEStepComponent*>	m_MC;	//!< array of model components active during this step

	DECLARE_PARAM_LIST();
};
