/*****************************************************************//**
 * \file   FEModel.h
 * \brief  This class define the FE model.
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once
#include "femcore/FEObjectBase.h"
#include "femcore/FEMesh.h"
#include "femcore/FELinearConstraintManager.h"
#include "femcore/FENLConstraint.h"
#include "femcore/FEGlobalData.h"
#include "basicio/DataStore.h"
#include "basicio/FEPlotDataStore.h"

#include <memory>

 // helper class for managing global (user-defined) variables.
class FEGlobalVariable
{
public:
    double v;
    std::string name;
};

/**
 * 定义整个求解模型. 包含所有的FEM组件，此类可进一步精简，定义每个组件的Manager
 */
class FEM_EXPORT FEModel: public FEObjectBase
{
public:
	FEModel(void);
	virtual ~FEModel(void);

	virtual bool Init() override;
	virtual bool Solve();
	bool IsSolved() const;

public:
	
	FEMesh& GetMesh();
	FELinearConstraintManager& GetLinearConstraintManager();

	bool InitBCs();
	bool InitMesh();

	//! Build the matrix profile for this model
	/**
	 * @~English
	 * @brief brief-description-about-BuildMatrixProfile .
	 * @param[??] G brief-description-about-G .
	 * @param[??] breset brief-description-about-breset .
	 * @return void brief-description-about-void .
	 */
	virtual void BuildMatrixProfile(FEGlobalMatrix& G, bool breset);

public:	

	void AddLoadController(FELoadController* plc);
	void ReplaceLoadController(int n, FELoadController* plc);
	FELoadController* GetLoadController(int i);

	int LoadControllers() const;

	void AttachLoadController(FEParam* p, int lc);
	void AttachLoadController(FEParam* p, FELoadController* plc);

	bool DetachLoadController(FEParam* p);
	FELoadController* GetLoadController(FEParam* p);

	//--- Material functions ---
public: 

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm);

	//! get the number of materials
	int Materials();

	//! return a pointer to a material
	FEMaterial* GetMaterial(int i);

	//! find a material based on its index
	FEMaterial* FindMaterial(int nid);

	//! find a material based on its name
	FEMaterial* FindMaterial(const std::string& matName);

	//! material initialization
	bool InitMaterials();

	//! material validation
	bool ValidateMaterials();

public:
	// Boundary conditions
	int BoundaryConditions() const;
	FEBoundaryCondition* BoundaryCondition(int i);
	void AddBoundaryCondition(FEBoundaryCondition* bc);
	void ClearBoundaryConditions();

	// initial conditions
	int InitialConditions();
	FEInitialCondition* InitialCondition(int i);
	void AddInitialCondition(FEInitialCondition* pbc);

public: // --- Analysis steps functions ---

	//! retrieve the number of steps
	int Steps();

	//! clear the steps
	void ClearSteps();

	//! Add an analysis step
	void AddStep(FEAnalysis* pstep);

	//! Get a particular step
	FEAnalysis* GetStep(int i);

	//! Get the current step
	FEAnalysis* GetCurrentStep();
	const FEAnalysis* GetCurrentStep() const;

	//! Set the current step index
	int GetCurrentStepIndex() const;

	//! Set the current step
	void SetCurrentStep(FEAnalysis* pstep);

	//! Set the current step index
	void SetCurrentStepIndex(int n);

	//! Get the current time
	FETimeInfo& GetTime();
	
	//! Get the start time
	double GetStartTime() const;

	//! Set the start time
	void SetStartTime(double t);

	//! Get the current time
	double GetCurrentTime() const;

	
	//! Set the current time
	void SetCurrentTime(double t);

	//! set the current time step
	void SetCurrentTimeStep(double dt);

public: // --- Contact interface functions ---

	//! return number of surface pair constraints
	int SurfacePairConstraints();

	//! retrive a surface pair interaction
	FESurfacePairConstraint* SurfacePairConstraint(int i);

	//! Add a surface pair constraint
	void AddSurfacePairConstraint(FESurfacePairConstraint* pci);

	//! Initializes contact data
	bool InitContact();

public: // --- Nonlinear constraints functions ---

	//! return number of nonlinear constraints
	int NonlinearConstraints();

	//! retrieve a nonlinear constraint
	FENLConstraint* NonlinearConstraint(int i);

	//! add a nonlinear constraint
	void AddNonlinearConstraint(FENLConstraint* pnlc);

	//! Initialize constraint data
	bool InitConstraints();

public:	// --- Model Loads ----
	//! return the number of model loads
	int ModelLoads();

	//! retrieve a model load
	FEModelLoad* ModelLoad(int i);

	//! Add a model load
	void AddModelLoad(FEModelLoad* pml);

	//! initialize model loads
	bool InitModelLoads();

public: // --- parameter functions ---

	//! evaluate all load controllers at some time
	void EvaluateLoadControllers(double time);

	// evaluate all mesh data
	void EvaluateDataGenerators(double time);

	//! evaluate all load parameters
	virtual bool EvaluateLoadParameters();

	//! return a reference to the named parameter 
public:	// --- Miscellaneous routines ---

	//! call the callback function
	//! This function returns fals if the run is to be aborted
	bool DoCallback(unsigned int nevent);

	//! I'd like to place the list of DOFS inside the model.
	//! As a first step, all classes that have access to the model
	//! should get the DOFS from this function
	DOFS& GetDOFS();

	//! Get the index of a DOF
	int GetDOFIndex(const char* sz) const;
	int GetDOFIndex(const char* szvar, int n) const;

	//! serialize data for restarts
	void Serialize(DumpStream& ar) override;

	//! This is called to serialize geometry.
	//! Derived classes can override this
	virtual void SerializeGeometry(DumpStream& ar);

	//! set the active module
	void SetActiveModule(const std::string& moduleName);

	//! get the module name
	string GetModuleName() const;

public: // Global data
	void AddGlobalData(FEGlobalData* psd);
	FEGlobalData* GetGlobalData(int i);
	FEGlobalData* FindGlobalData(const char* szname);
	int FindGlobalDataIndex(const char* szname);
	int GlobalDataItems();

	// get/set global data
	void SetGlobalConstant(const string& s, double v);
	double GetGlobalConstant(const string& s);

	int GlobalVariables() const;
	void AddGlobalVariable(const std::string& s, double v);
	const FEGlobalVariable& GetGlobalVariable(int n);

public: // Data retrieval

	// get nodal dof data
	bool GetNodeData(int dof, vector<double>& data);

	//! return the data store
	DataStore& GetDataStore();

	//! return plot data
	FEPlotDataStore& GetPlotDataStore();
	const FEPlotDataStore& GetPlotDataStore() const;

public:
	// reset all the timers
	void ResetAllTimers();

	// return total number of timers
	int Timers();

	// return a timer by index
	Timer* GetTimer(int i);

	// get the number of calls to Update()
	int UpdateCounter() const;

	// this can be used to change the update counter
	void IncrementUpdateCounter();

public:
	void SetUnits(const char* szunits);
	const char* GetUnits() const;

private:
	class Impl;
	std::unique_ptr<Impl> m_imp;
	DECLARE_PARAM_LIST();
};
