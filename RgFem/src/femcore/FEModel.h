/*****************************************************************//**
 * \file   FEModel.h
 * \brief  This class define the FE model.
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"

#include <memory>

/**
 * 定义整个求解模型.整个
 */
class FE_EXPORT FEModel : public FEObjectBase
{
public:
	FEModel(void);
	virtual ~FEModel(void);

	// Initialization
	virtual bool Init() override;

	// solve the model
	virtual bool Solve();

	//! Call this function whenever the geometry of the model has changed.
	virtual void Update();

	//! will return true if the model solved succussfully
	bool IsSolved() const;

public:
	// get the FE mesh
	FEMesh& GetMesh();

	// get the linear constraint manager
	FELinearConstraintManager& GetLinearConstraintManager();

	//! Validate BC's
	bool InitBCs();

	//! Initialize the mesh
	bool InitMesh();

	//! Initialize shells
	virtual void InitShells();

	//! Build the matrix profile for this model
	/**
	 * @~English
	 * @brief brief-description-about-BuildMatrixProfile .
	 * @param[??] G brief-description-about-G .
	 * @param[??] breset brief-description-about-breset .
	 * @return void brief-description-about-void .
	 */
	virtual void BuildMatrixProfile(FEGlobalMatrix& G, bool breset);

public:	// --- Load controller functions ----

	//! Add a load controller to the model
	void AddLoadController(FELoadController* plc);

	//! replace a load controller
	void ReplaceLoadController(int n, FELoadController* plc);

	//! get a load controller
	FELoadController* GetLoadController(int i);

	//! get the number of load controllers
	int LoadControllers() const;

	//! Attach a load controller to a parameter
	void AttachLoadController(FEParam* p, int lc);
	void AttachLoadController(FEParam* p, FELoadController* plc);

	//! Detach a load controller from a parameter
	bool DetachLoadController(FEParam* p);

	//! Get a load controller for a parameter (returns null if the param is not under load control)
	FELoadController* GetLoadController(FEParam* p);

public:	// --- mesh data generators ---

	//! Add a mesh data generator to the model
	void AddMeshDataGenerator(FEMeshDataGenerator* pmd);

	//! get a mesh data generator
	FEMeshDataGenerator* GetMeshDataGenerator(int i);

	//! get the number of mesh data generators
	int MeshDataGenerators() const;

public: // --- Material functions ---

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

public:	// --- Mesh adaptors ---
	//! return number of mesh adaptors
	int MeshAdaptors();

	//! retrieve a mesh adaptors
	FEMeshAdaptor* MeshAdaptor(int i);

	//! add a mesh adaptor
	void AddMeshAdaptor(FEMeshAdaptor* meshAdaptor);

public: // --- parameter functions ---

	//! evaluate all load controllers at some time
	void EvaluateLoadControllers(double time);

	// evaluate all mesh data
	void EvaluateDataGenerators(double time);

	//! evaluate all load parameters
	virtual bool EvaluateLoadParameters();

	//! Find a model parameter
	FEParam* FindParameter(const ParamString& s) override;

	//! return a reference to the named parameter
	virtual FEParamValue GetParameterValue(const ParamString& param);

	//! Find property 
	//! Note: Can't call this FindProperty, since this is already defined in base class
	FECoreBase* FindComponent(const ParamString& prop);

	//! Set the print parameters flag
	void SetPrintParametersFlag(bool b);

	//! Get the print parameter flag
	bool GetPrintParametersFlag() const;

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

public:
	//! Log a message
	void Logf(int ntag, const char* msg, ...);
	void BlockLog();
	void UnBlockLog();
	bool LogBlocked() const;

public:
	// Derived classes can use this to implement the actual logging mechanism
	virtual void Log(int ntag, const char* msg);

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
	void AddGlobalVariable(const string& s, double v);
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

protected:
	FEParamValue GetMeshParameter(const ParamString& paramString);

private:
	class Impl;
	std::unique_ptr<Impl> m_imp;
	DECLARE_FECORE_CLASS();
};

FECORE_API FECoreBase* CopyFEBioClass(FECoreBase* pc, FEModel* fem);
