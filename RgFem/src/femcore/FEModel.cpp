#include "FEModel.h"

#include "basicio/DataStore.h"
#include "basicio/DumpMemStream.h"
#include "basicio/DumpStream.h"
#include "basicio/FEPlotDataStore.h"
#include "Domain/FEDomain.h"
#include "FEBoundaryCondition.h"
#include "FEException.h"
#include "FEGlobalData.h"
#include "FEInitialCondition.h"
#include "FELinearConstraintManager.h"
#include "FELoadController.h"
#include "femcore/Callback.h"
#include "femcore/FEAnalysis/FEAnalysis.h"
#include "femcore/FEMesh.h"
#include "FEModelLoad.h"
#include "FEModelParam.h"
#include "FEModule.h"
#include "FENLConstraint.h"
#include "FESolidModule.h"
#include "FESurfaceLoad.h"
#include "FESurfacePairConstraint.h"
#include "logger/log.h"
#include "materials/FEMaterial.h"
#include "Timer.h"

#include <map>
#include <stdarg.h>
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
// Impl class for the FEModel class
class FEModel::Impl
{
public:
    struct LoadParam
    {
        FEParam* param;
        int lc;

        double m_scl;
        Vector3d m_vscl;

        LoadParam()
        {
            m_scl = 1.0;
            m_vscl = Vector3d(0, 0, 0);
        }

        void Serialize(DumpStream& ar)
        {
            ar& lc;
            ar& m_scl& m_vscl;

            if (ar.IsShallow() == false)
            {
                // we can't save the FEParam* directly, so we need to store meta data and try to find it on loading
                if (ar.IsSaving())
                {
                    FEObjectBase* pc = dynamic_cast<FEObjectBase*>(param->owner());
                    assert(pc);
                    ar << pc;
                    ar << param->name();
                }
                else
                {
                    FEObjectBase* pc = nullptr;
                    ar >> pc;
                    assert(pc);

                    char name[256] = {0};
                    ar >> name;

                    param = pc->FindParameter(name);
                    assert(param);
                }
            }
            else
                param = nullptr;
        }
    };

public:
    Impl(FEModel* fem)
        : m_fem(fem)
        , m_mesh(fem)
        , m_dmp(*fem)
    {
        // --- Analysis Data ---
        mp_CurStep = 0;
        m_nStep = -1;
        m_ftime0 = 0;

        m_nupdates = 0;
        m_bsolved = false;
        m_block_log = false;
        m_printParams = false;
        m_meshUpdate = false;

        // create the linear constraint manager
        m_LCM = new FELinearConstraintManager(fem);

        m_timers.resize(7);
    }

public:                     // TODO: Find a better place for these parameters
    FETimeInfo m_timeInfo;  //!< current time value
    double m_ftime0;        //!< start time of current step

    bool m_block_log;
    int m_nupdates;  //!< number of calls to FEModel::Update

public:
    std::vector<FEMaterial*> m_MAT;              //!< array of materials
    std::vector<FEBoundaryCondition*> m_BC;      //!< boundary conditions
    std::vector<FEModelLoad*> m_ML;              //!< model loads
    std::vector<FEInitialCondition*> m_IC;       //!< initial conditions
    std::vector<FESurfacePairConstraint*> m_CI;  //!< contact interface array
    std::vector<FENLConstraint*> m_NLC;          //!< nonlinear constraints
    std::vector<FELoadController*> m_LC;         //!< load controller data
    std::vector<FEAnalysis*> m_Step;             //!< array of analysis steps
    std::vector<LoadParam> m_Param;              //!< list of parameters controller by load controllers
    std::vector<Timer> m_timers;                 // list of timers

public:
    FEAnalysis* mp_CurStep;  //!< pointer to current analysis step
    int m_nStep;             //!< current index of analysis step
    bool m_printParams;      //!< print parameters
    bool m_meshUpdate;       //!< mesh update flag
    std::string m_units;     // units string

public:
    FEModel* m_fem;
    std::string m_moduleName;
    bool m_bsolved;              // solved flag
    DOFS m_dofs;                 // meta info of degree of freedoms in this model
    FEMesh m_mesh;               // the one and only FE mesh
    FELinearConstraintManager* m_LCM;
    DataStore m_dataStore;       //!< the data store used for data logging
    FEPlotDataStore m_plotData;  //!< Output request for plot file
    DumpMemStream m_dmp;         // only used by incremental solver
};

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEModel, FEObjectBase)
// model parameters
ADD_PARAMETER(m_imp->m_timeInfo.currentTime, "time");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEModel::FEModel(void)
    : FEObjectBase()
    , m_imp(new FEModel::Impl(this))
{
    // set the name
    SetName("fem");

    // reset all timers
    ResetAllTimers();
}

//-----------------------------------------------------------------------------
//! Delete all dynamically allocated data
FEModel::~FEModel(void)
{
    Clear();
}

//-----------------------------------------------------------------------------
//! return the data store
DataStore& FEModel::GetDataStore()
{
    return m_imp->m_dataStore;
}

//-----------------------------------------------------------------------------
FEPlotDataStore& FEModel::GetPlotDataStore()
{
    return m_imp->m_plotData;
}

//-----------------------------------------------------------------------------
const FEPlotDataStore& FEModel::GetPlotDataStore() const
{
    return m_imp->m_plotData;
}

//-----------------------------------------------------------------------------
//! will return true if the model solved succussfully
bool FEModel::IsSolved() const
{
    return m_imp->m_bsolved;
}

//-----------------------------------------------------------------------------
void FEModel::Clear()
{
    // clear dofs
    m_imp->m_dofs.Reset();

    // clear all properties
    for (FEMaterial* mat : m_imp->m_MAT)
        delete mat;
    m_imp->m_MAT.clear();
    for (FEBoundaryCondition* bc : m_imp->m_BC)
        delete bc;
    m_imp->m_BC.clear();
    for (FEModelLoad* ml : m_imp->m_ML)
        delete ml;
    m_imp->m_ML.clear();
    for (FEInitialCondition* ic : m_imp->m_IC)
        delete ic;
    m_imp->m_IC.clear();
    for (FESurfacePairConstraint* ci : m_imp->m_CI)
        delete ci;
    m_imp->m_CI.clear();
    for (FENLConstraint* nlc : m_imp->m_NLC)
        delete nlc;
    m_imp->m_NLC.clear();
    for (FELoadController* lc : m_imp->m_LC)
        delete lc;
    m_imp->m_LC.clear();
    /*for (FEMeshDataGenerator* md : m_imp->m_MD)
        delete md;
    m_imp->m_MD.clear();*/
    for (FEAnalysis* step : m_imp->m_Step)
        delete step;
    m_imp->m_Step.clear();

    // global data
    /* for (size_t i = 0; i < m_imp->m_GD.size(); ++i)
         delete m_imp->m_GD[i];
     m_imp->m_GD.clear();
     m_imp->m_Const.clear();*/

    // global variables (TODO: Should I delete the corresponding parameters?)
    /*for (size_t i = 0; i < m_imp->m_Var.size(); ++i)
        delete m_imp->m_Var[i];
    m_imp->m_Var.clear();*/

    // clear the linear constraints
    if (m_imp->m_LCM)
        m_imp->m_LCM->Clear();

    // clear the mesh
    m_imp->m_mesh.Clear();

    // clear load parameters
    /*m_imp->m_Param.clear();*/
}

//-----------------------------------------------------------------------------
//! set the module name
void FEModel::SetActiveModule(const std::string& moduleName)
{
    m_imp->m_moduleName = moduleName;
    // FECoreKernel& fecore = FECoreKernel::GetInstance();
    // fecore.SetActiveModule(moduleName.c_str());
    // FEModule* pmod = fecore.GetActiveModule();
    // pmod->InitModel(this);

    FEModule* pmod = new FESolidModule();
    pmod->InitModel(this);
}

//-----------------------------------------------------------------------------
//! get the module name
string FEModel::GetModuleName() const
{
    return m_imp->m_moduleName;
}

//-----------------------------------------------------------------------------
int FEModel::BoundaryConditions() const
{
    return (int)m_imp->m_BC.size();
}

//-----------------------------------------------------------------------------
FEBoundaryCondition* FEModel::BoundaryCondition(int i)
{
    return m_imp->m_BC[i];
}

//-----------------------------------------------------------------------------
void FEModel::AddBoundaryCondition(FEBoundaryCondition* pbc)
{
    m_imp->m_BC.push_back(pbc);
}

//-----------------------------------------------------------------------------
void FEModel::ClearBoundaryConditions()
{
    m_imp->m_BC.clear();
}

//-----------------------------------------------------------------------------
int FEModel::InitialConditions()
{
    return (int)m_imp->m_IC.size();
}

//-----------------------------------------------------------------------------
FEInitialCondition* FEModel::InitialCondition(int i)
{
    return m_imp->m_IC[i];
}

//-----------------------------------------------------------------------------
void FEModel::AddInitialCondition(FEInitialCondition* pbc)
{
    m_imp->m_IC.push_back(pbc);
}

//-----------------------------------------------------------------------------
//! retrieve the number of steps
int FEModel::Steps() const
{
    return (int)m_imp->m_Step.size();
}

//-----------------------------------------------------------------------------
//! clear the steps
void FEModel::ClearSteps()
{
    m_imp->m_Step.clear();
}

//-----------------------------------------------------------------------------
//! Add an analysis step
void FEModel::AddStep(FEAnalysis* pstep)
{
    m_imp->m_Step.push_back(pstep);
}

//-----------------------------------------------------------------------------
//! Get a particular step
FEAnalysis* FEModel::GetStep(int i)
{
    return m_imp->m_Step[i];
}

//-----------------------------------------------------------------------------
//! Get the current step
FEAnalysis* FEModel::GetCurrentStep()
{
    return m_imp->mp_CurStep;
}
const FEAnalysis* FEModel::GetCurrentStep() const
{
    return m_imp->mp_CurStep;
}

//-----------------------------------------------------------------------------
//! Set the current step
void FEModel::SetCurrentStep(FEAnalysis* pstep)
{
    m_imp->mp_CurStep = pstep;
}

//-----------------------------------------------------------------------------
//! Set the current step index
int FEModel::GetCurrentStepIndex() const
{
    return m_imp->m_nStep;
}

//-----------------------------------------------------------------------------
//! Set the current step index
void FEModel::SetCurrentStepIndex(int n)
{
    m_imp->m_nStep = n;
}

//-----------------------------------------------------------------------------
//! return number of surface pair interactions
int FEModel::SurfacePairConstraints()
{
    return (int)m_imp->m_CI.size();
}

//-----------------------------------------------------------------------------
//! retrive a surface pair interaction
FESurfacePairConstraint* FEModel::SurfacePairConstraint(int i)
{
    return m_imp->m_CI[i];
}

//-----------------------------------------------------------------------------
//! Add a surface pair interaction
void FEModel::AddSurfacePairConstraint(FESurfacePairConstraint* pci)
{
    m_imp->m_CI.push_back(pci);
}

//-----------------------------------------------------------------------------
//! return number of nonlinear constraints
int FEModel::NonlinearConstraints()
{
    return (int)m_imp->m_NLC.size();
}

//-----------------------------------------------------------------------------
//! retrieve a nonlinear constraint
FENLConstraint* FEModel::NonlinearConstraint(int i)
{
    return m_imp->m_NLC[i];
}

//-----------------------------------------------------------------------------
//! add a nonlinear constraint
void FEModel::AddNonlinearConstraint(FENLConstraint* pnlc)
{
    m_imp->m_NLC.push_back(pnlc);
}

//-----------------------------------------------------------------------------
//! return the number of model loads
int FEModel::ModelLoads()
{
    return (int)m_imp->m_ML.size();
}

//-----------------------------------------------------------------------------
//! retrieve a model load
FEModelLoad* FEModel::ModelLoad(int i)
{
    return m_imp->m_ML[i];
}

//-----------------------------------------------------------------------------
//! Add a model load
void FEModel::AddModelLoad(FEModelLoad* pml)
{
    m_imp->m_ML.push_back(pml);
}

//-----------------------------------------------------------------------------
// get the FE mesh
FEMesh& FEModel::GetMesh()
{
    return m_imp->m_mesh;
}

//-----------------------------------------------------------------------------
FELinearConstraintManager& FEModel::GetLinearConstraintManager()
{
    return *m_imp->m_LCM;
}

//-----------------------------------------------------------------------------
bool FEModel::Init()
{
    // make sure there is something to do
    if (m_imp->m_Step.size() == 0)
        return false;

    // intitialize time
    FETimeInfo& tp = GetTime();
    tp.currentTime = 0;
    m_imp->m_ftime0 = 0;

    // evaluate all load controllers at the initial time
    for (int i = 0; i < LoadControllers(); ++i)
    {
        FELoadController* plc = m_imp->m_LC[i];
        if (plc->Init() == false)
        {
            std::string s = plc->GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Load controller %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
        plc->Evaluate(0);
    }

    // evaluate all mesh data generators
    /*  for (int i = 0; i < MeshDataGenerators(); ++i)
      {
          FEMeshDataGenerator* pmd = m_imp->m_MD[i];
          if (pmd->Init() == false)
          {
              std::string s = pmd->GetName();
              const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
              feLogError("Node data generator %d (%s) failed to initialize", i + 1, sz);
              return false;
          }
          pmd->Evaluate(0);
      }*/

    // check step data
    for (int i = 0; i < (int)m_imp->m_Step.size(); ++i)
    {
        FEAnalysis& step = *m_imp->m_Step[i];
        if (step.Init() == false)
        {
            std::string s = step.GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Step %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }

    // validate BC's
    if (InitBCs() == false)
        return false;

    // initialize material data
    // NOTE: This must be called after the rigid system is initialiazed since the rigid materials will
    //       reference the rigid bodies
    if (InitMaterials() == false)
        return false;

    // initialize mesh data
    // NOTE: this must be done AFTER the elements have been assigned material point data !
    // this is because the mesh data is reset
    // TODO: perhaps I should not reset the mesh data during the initialization
    if (InitMesh() == false)
        return false;

    // initialize model loads
    // NOTE: This must be called after the InitMaterials since the COM of the rigid bodies
    //       are set in that function.
    if (InitModelLoads() == false)
        return false;

    // initialize contact data
    if (InitContact() == false)
        return false;

    // initialize nonlinear constraints
    if (InitConstraints() == false)
        return false;

    // evaluate all load parameters
    // Do this last in case any model components redefined their load curves.
    /*if (EvaluateLoadParameters() == false)
        return false;*/

    // activate all permanent dofs
    // Activate();
    for (auto pbc : m_imp->m_BC)
    {
        if (pbc)
        {
            pbc->Activate();
        }
    }

    bool ret = false;
    try
    {
        ret = DoCallback(CB_INIT);
    }
    catch (std::exception c)
    {
        ret = false;
        feLogError(c.what());
    }

    // do the callback
    return ret;
}

//-----------------------------------------------------------------------------
// get the number of calls to Update()
int FEModel::UpdateCounter() const
{
    return m_imp->m_nupdates;
}

//-----------------------------------------------------------------------------
void FEModel::IncrementUpdateCounter()
{
    m_imp->m_nupdates++;
}

//-----------------------------------------------------------------------------
void FEModel::Update()
{
    TRACK_TIME(TimerID::Timer_Update);

    // update model counter
    m_imp->m_nupdates++;

    // update mesh
    FEMesh& mesh = GetMesh();
    const FETimeInfo& tp = GetTime();
    mesh.Update(tp);

    // set the mesh update flag to false
    // If any load sets this to true, the
    // mesh will also be update after the loads are updated
    m_imp->m_meshUpdate = false;

    int nvel = BoundaryConditions();
    for (int i = 0; i < nvel; ++i)
    {
        FEBoundaryCondition& bc = *BoundaryCondition(i);
        if (bc.IsActive())
            bc.UpdateModel();
    }

    // update all model loads
    for (int i = 0; i < ModelLoads(); ++i)
    {
        FEModelLoad* pml = ModelLoad(i);
        if (pml && pml->IsActive())
            pml->Update();
    }

    // update all paired-interfaces
    for (int i = 0; i < SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* psc = SurfacePairConstraint(i);
        if (psc && psc->IsActive())
            psc->Update();
    }

    // update all constraints
    for (int i = 0; i < NonlinearConstraints(); ++i)
    {
        FENLConstraint* pc = NonlinearConstraint(i);
        if (pc && pc->IsActive())
            pc->Update();
    }

    // some of the loads may alter the prescribed dofs, so we update the mesh again
    if (m_imp->m_meshUpdate)
    {
        mesh.Update(tp);
        m_imp->m_meshUpdate = false;
    }

    // do the callback
    DoCallback(CB_MODEL_UPDATE);
}

//-----------------------------------------------------------------------------
//! See if the BC's are setup correctly.
bool FEModel::InitBCs()
{
    // check the IC's
    int NIC = InitialConditions();
    for (int i = 0; i < NIC; ++i)
    {
        FEInitialCondition* pic = InitialCondition(i);
        if (pic->Init() == false)
        {
            std::string s = pic->GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Initial condition %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }

    // check the BC's
    int NBC = BoundaryConditions();
    for (int i = 0; i < NBC; ++i)
    {
        FEBoundaryCondition* pbc = BoundaryCondition(i);
        if (pbc->Init() == false)
        {
            std::string s = pbc->GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Boundary condition %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }

    // check the linear constraints
    if (GetLinearConstraintManager().Initialize() == false)
        return false;

    return true;
}

//-----------------------------------------------------------------------------
void FEModel::AddMaterial(FEMaterial* pm)
{
    m_imp->m_MAT.push_back(pm);
}

//-----------------------------------------------------------------------------
//! get the number of materials
int FEModel::Materials()
{
    return (int)m_imp->m_MAT.size();
}

//-----------------------------------------------------------------------------
//! return a pointer to a material
FEMaterial* FEModel::GetMaterial(int i)
{
    return m_imp->m_MAT[i];
}

//-----------------------------------------------------------------------------
FEMaterial* FEModel::FindMaterial(int nid)
{
    for (int i = 0; i < Materials(); ++i)
    {
        FEMaterial* pm = GetMaterial(i);
        if (pm->GetID() == nid)
            return pm;
    }
    return 0;
}

//-----------------------------------------------------------------------------
FEMaterial* FEModel::FindMaterial(const std::string& matName)
{
    for (int i = 0; i < Materials(); ++i)
    {
        FEMaterial* mat = GetMaterial(i);
        if (mat->GetName() == matName)
            return mat;
    }
    return 0;
}

//-----------------------------------------------------------------------------
//! Initialize material data (This also does an initial validation).
bool FEModel::InitMaterials()
{
    // initialize material data
    for (int i = 0; i < Materials(); ++i)
    {
        // get the material
        FEMaterial* pmat = GetMaterial(i);

        // initialize material data
        if (pmat->Init() == false)
        {
            feLogError("Failed initializing material %d (name=\"%s\")", i + 1, pmat->GetName().c_str());
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! validate material data
bool FEModel::ValidateMaterials()
{
    // initialize material data
    for (int i = 0; i < Materials(); ++i)
    {
        // get the material
        FEMaterial* pmat = GetMaterial(i);

        // initialize material data
        if (pmat->Validate() == false)
        {
            feLogError("Failed validating material %d (name=\"%s\")", i + 1, pmat->GetName().c_str());
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! Add a loadcurve to the model
void FEModel::AddLoadController(FELoadController* plc)
{
    m_imp->m_LC.push_back(plc);
}

//-----------------------------------------------------------------------------
void FEModel::ReplaceLoadController(int n, FELoadController* plc)
{
    assert((n >= 0) && (n < LoadControllers()));
    delete m_imp->m_LC[n];
    m_imp->m_LC[n] = plc;
}

//-----------------------------------------------------------------------------
//! get a loadcurve
FELoadController* FEModel::GetLoadController(int i)
{
    return m_imp->m_LC[i];
}

//-----------------------------------------------------------------------------
//! get the number of loadcurves
int FEModel::LoadControllers() const
{
    return (int)m_imp->m_LC.size();
}

//-----------------------------------------------------------------------------
//! Attach a load controller to a parameter
void FEModel::AttachLoadController(FEParam* param, int lc)
{
    Impl::LoadParam lp;
    lp.param = param;
    lp.lc = lc;

    switch (param->type())
    {
        case FE_PARAM_DOUBLE:
            lp.m_scl = param->value<double>();
            break;
        case FE_PARAM_VEC3D:
            lp.m_vscl = param->value<Vector3d>();
            break;
    }

    m_imp->m_Param.push_back(lp);
}

//-----------------------------------------------------------------------------
void FEModel::AttachLoadController(FEParam* p, FELoadController* plc)
{
    AttachLoadController(p, plc->GetID());
}

//-----------------------------------------------------------------------------
//! Detach a load controller from a parameter
bool FEModel::DetachLoadController(FEParam* p)
{
    for (int i = 0; i < (int)m_imp->m_Param.size(); ++i)
    {
        Impl::LoadParam& pi = m_imp->m_Param[i];
        if (pi.param == p)
        {
            m_imp->m_Param.erase(m_imp->m_Param.begin() + i);
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
//! Get a load controller for a parameter (returns null if the param is not under load control)
FELoadController* FEModel::GetLoadController(FEParam* p)
{
    /*for (int i = 0; i < (int)m_imp->m_Param.size(); ++i)
    {
        Impl::LoadParam& pi = m_imp->m_Param[i];
        if (pi.param == p)
        {
            return (pi.lc >= 0 ? GetLoadController(pi.lc) : nullptr);
        }
    }*/
    return nullptr;
}

//-----------------------------------------------------------------------------
//! Initialize rigid force data
bool FEModel::InitModelLoads()
{
    // call the Init() function of all rigid forces
    for (int i = 0; i < ModelLoads(); ++i)
    {
        FEModelLoad& FC = *ModelLoad(i);
        if (FC.Init() == false)
        {
            std::string s = FC.GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Load %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Does one-time initialization of the Mesh data. Call FEMesh::Reset for resetting
//! the mesh data.
bool FEModel::InitMesh()
{
    FEMesh& mesh = GetMesh();

    // find and remove isolated vertices
    int ni = mesh.RemoveIsolatedVertices();
    if (ni != 0)
    {
        if (ni == 1)
            feLogWarning("%d isolated vertex removed.", ni);
        else
            feLogWarning("%d isolated vertices removed.", ni);
    }


    //// Initialize shell data
    //// This has to be done before the domains are initialized below
    // InitShells();

    // reset data
    // TODO: Not sure why this is here
    try
    {
        mesh.Reset();
    }
    catch (NegativeJacobian e)
    {
        feLogError("Negative jacobian detected during mesh initialization.");
        return false;
    }

    // initialize all domains
    // Initialize shell domains first (in order to establish SSI)
    // TODO: I'd like to move the initialization of the SSI to InitShells, but I can't
    //       do that because FESSIShellDomain::FindSSI depends on the FEDomain::m_Node array which is
    //       initialized in FEDomain::Init.
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.Class() == FE_DOMAIN_SHELL)
            if (dom.Init() == false)
                return false;
    }

    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.Class() != FE_DOMAIN_SHELL)
            if (dom.Init() == false)
                return false;
    }

    //// initialize surfaces
    // for (int i = 0; i < mesh.Surfaces(); ++i)
    //{
    //     if (mesh.Surface(i).Init() == false)
    //         return false;
    // }

    // All done
    return true;
}

//-----------------------------------------------------------------------------
//! Initializes contact data
bool FEModel::InitContact()
{
    // loop over all contact interfaces
    for (int i = 0; i < SurfacePairConstraints(); ++i)
    {
        // get the contact interface
        FESurfacePairConstraint& ci = *SurfacePairConstraint(i);

        // initializes contact interface data
        if (ci.Init() == false)
        {
            std::string s = ci.GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Contact %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! Initialize the nonlinear constraints.
//! This function is called during model initialization (\sa FEModel::Init)
bool FEModel::InitConstraints()
{
    for (int i = 0; i < NonlinearConstraints(); ++i)
    {
        FENLConstraint* plc = NonlinearConstraint(i);

        // initialize
        if (plc->Init() == false)
        {
            std::string s = plc->GetName();
            const char* sz = (s.empty() ? "<unnamed>" : s.c_str());
            feLogError("Nonlinear constraint %d (%s) failed to initialize", i + 1, sz);
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
//! This function solves the FE problem by calling the solve method for each step.
bool FEModel::Solve()
{
    TRACK_TIME(Timer_ModelSolve);

    // error flag
    bool bOk = true;

    // loop over all analysis steps
    // Note that we don't necessarily from step 0 as user can use restart~
    for (size_t iStep = m_imp->m_nStep; iStep < Steps(); ++iStep)
    {
        // set the current analysis step
        m_imp->m_nStep = iStep;
        m_imp->mp_CurStep = m_imp->m_Step[(int)iStep];

        // intitialize step data
        if (m_imp->mp_CurStep->Activate() == false)
        {
            bOk = false;
            break;
        }

        DoCallback(CB_STEP_ACTIVE);

        // solve the analaysis step
        bOk = m_imp->mp_CurStep->Solve();

        if (iStep + 1 == Steps())  // the last step
        {
            m_imp->m_bsolved = bOk;
        }

        // do callbacks
        DoCallback(CB_STEP_SOLVED);

        // break if the step has failed
        if (!bOk)
            break;
    }

    // do the callbacks
    DoCallback(CB_SOLVED);

    return bOk;
}

//-----------------------------------------------------------------------------
//! Get the current time information.
FETimeInfo& FEModel::GetTime()
{
    return m_imp->m_timeInfo;
}

//-----------------------------------------------------------------------------
double FEModel::GetStartTime() const
{
    return m_imp->m_ftime0;
}

//-----------------------------------------------------------------------------
void FEModel::SetStartTime(double t)
{
    m_imp->m_ftime0 = t;
}

//-----------------------------------------------------------------------------
double FEModel::GetCurrentTime() const
{
    return m_imp->m_timeInfo.currentTime;
}

//-----------------------------------------------------------------------------
void FEModel::SetCurrentTime(double t)
{
    m_imp->m_timeInfo.currentTime = t;
}

//-----------------------------------------------------------------------------
void FEModel::SetCurrentTimeStep(double dt)
{
    FEAnalysis* step = GetCurrentStep();
    assert(step);
    if (step)
        step->m_dt = dt;
}

//-----------------------------------------------------------------------------
DOFS& FEModel::GetDOFS()
{
    return m_imp->m_dofs;
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* sz) const
{
    return m_imp->m_dofs.GetDOF(sz);
}

//-----------------------------------------------------------------------------
int FEModel::GetDOFIndex(const char* szvar, int n) const
{
    return m_imp->m_dofs.GetDOF(szvar, n);
}

//-----------------------------------------------------------------------------
// Call the callback function if there is one defined
//
bool FEModel::DoCallback(unsigned int nevent)
{
    try
    {
        // do the callbacks
        // Rango TODO:
        // bool bret = CallbackHandler::DoCallback(this, nevent);
        // return bret;
    }
    catch (ForceConversion)
    {
        throw;
    }
    catch (IterationFailure)
    {
        throw;
    }
    catch (DoRunningRestart)
    {
        throw;
    }
    catch (std::exception e)
    {
        throw;
    }
    catch (...)
    {
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
void FEModel::SetGlobalConstant(const string& s, double v)
{
    // m_imp->m_Const[s] = v;
    return;
}

//-----------------------------------------------------------------------------
double FEModel::GetGlobalConstant(const string& s)
{
    // return (m_imp->m_Const.count(s) ? m_imp->m_Const.find(s)->second : 0);
    return 0;
}

//-----------------------------------------------------------------------------
int FEModel::GlobalVariables() const
{
    // return (int)m_imp->m_Var.size();
    return 0;
}

//-----------------------------------------------------------------------------
void FEModel::AddGlobalVariable(const string& s, double v)
{
    FEGlobalVariable* var = new FEGlobalVariable;
    var->v = v;
    var->name = s;
    AddParameter(var->v, var->name.c_str());
    // m_imp->m_Var.push_back(var);
}

const FEGlobalVariable& FEModel::GetGlobalVariable(int n)
{
    // return *m_imp->m_Var[n];
    return FEGlobalVariable();
}

//-----------------------------------------------------------------------------
void FEModel::AddGlobalData(FEGlobalData* psd)
{
    // m_imp->m_GD.push_back(psd);
}

//-----------------------------------------------------------------------------
FEGlobalData* FEModel::GetGlobalData(int i)
{
    // return m_imp->m_GD[i];
    return nullptr;
}

//-----------------------------------------------------------------------------
FEGlobalData* FEModel::FindGlobalData(const char* szname)
{
    /*for (int i = 0; i < m_imp->m_GD.size(); ++i)
    {
        if (m_imp->m_GD[i]->GetName() == szname)
            return m_imp->m_GD[i];
    }*/
    return nullptr;
}

//-----------------------------------------------------------------------------
int FEModel::FindGlobalDataIndex(const char* szname)
{
    /*for (int i = 0; i < m_imp->m_GD.size(); ++i)
    {
        if (m_imp->m_GD[i]->GetName() == szname)
            return i;
    }*/
    return -1;
}

//-----------------------------------------------------------------------------
int FEModel::GlobalDataItems()
{
    // return (int)m_imp->m_GD.size();
    return 0;
}

//-----------------------------------------------------------------------------
// This function serializes data to a stream.
// This is used for running and cold restarts.
// void FEModel::Impl::Serialize(DumpStream& ar)
//{
//    if (ar.IsShallow())
//    {
//        // stream model data
//        ar& m_timeInfo;
//
//        // stream mesh
//        m_fem->SerializeGeometry(ar);
//
//        // serialize contact
//        ar& m_CI;
//
//        // serialize nonlinear constraints
//        ar& m_NLC;
//
//        // serialize step and solver data
//        ar& m_Step;
//    }
//    else
//    {
//        if (ar.IsLoading())
//            m_fem->Clear();
//
//        ar& m_moduleName;
//
//        if (ar.IsLoading())
//        {
//            FECoreKernel::GetInstance().SetActiveModule(m_moduleName.c_str());
//        }
//
//        ar& m_timeInfo;
//        ar& m_dofs;
//        ar& m_Const;
//        ar& m_GD;
//        ar& m_ftime0;
//        ar& m_bsolved;
//
//        // we have to stream materials before the mesh
//        ar& m_MAT;
//
//        // we have to stream the mesh before any boundary conditions
//        m_fem->SerializeGeometry(ar);
//
//        // stream all boundary conditions
//        ar& m_BC;
//        ar& m_ML;
//        ar& m_IC;
//        ar& m_CI;
//        ar& m_NLC;
//
//        // stream step data next
//        ar& m_nStep;
//        ar& m_Step;
//        ar& mp_CurStep;  // This must be streamed after m_Step
//
//        // serialize linear constraints
//        if (m_LCM)
//            m_LCM->Serialize(ar);
//
//        // serialize data generators
//        ar& m_MD;
//
//        // load controllers and load parameters are streamed last
//        // since they can depend on other model parameters.
//        ar& m_LC;
//        ar& m_Param;
//    }
//}

//-----------------------------------------------------------------------------
//! This is called to serialize geometry.
//! Derived classes can override this
void FEModel::SerializeGeometry(DumpStream& ar)
{
    ar & m_imp->m_mesh;
}

//-----------------------------------------------------------------------------
// This function serializes data to a stream.
// This is used for running and cold restarts.
void FEModel::Serialize(DumpStream& ar)
{
    TRACK_TIME(TimerID::Timer_Update);

    // m_imp->Serialize(ar);
    DoCallback(ar.IsSaving() ? CB_SERIALIZE_SAVE : CB_SERIALIZE_LOAD);
}

//-----------------------------------------------------------------------------
void FEModel::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
    FEAnalysis* pstep = GetCurrentStep();
    FESolver* solver = pstep->GetFESolver();
    solver->BuildMatrixProfile(G, breset);
}

//-----------------------------------------------------------------------------
bool FEModel::GetNodeData(int ndof, vector<double>& data)
{
    // get the dofs
    DOFS& dofs = GetDOFS();

    // make sure the dof index is valid
    if ((ndof < 0) || (ndof >= dofs.GetTotalDOFS()))
        return false;

    // get the mesh and number of nodes
    FEMesh& mesh = GetMesh();
    int N = mesh.Nodes();

    // make sure data is correct size
    data.resize(N, 0.0);

    // loop over all nodes
    for (int i = 0; i < N; ++i)
    {
        FENode& node = mesh.Node(i);
        data[i] = node.get(ndof);
    }

    return true;
}

//-----------------------------------------------------------------------------
// reset all the timers
void FEModel::ResetAllTimers()
{
    for (size_t i = 0; i < m_imp->m_timers.size(); ++i)
    {
        Timer& ti = m_imp->m_timers[i];
        ti.reset();
    }
}

//-----------------------------------------------------------------------------
int FEModel::Timers()
{
    return (int)m_imp->m_timers.size();
}

//-----------------------------------------------------------------------------
Timer* FEModel::GetTimer(int i)
{
    return &(m_imp->m_timers[i]);
}

//-----------------------------------------------------------------------------
void FEModel::SetUnits(const char* szunits)
{
    if (szunits)
        m_imp->m_units = szunits;
    else
        m_imp->m_units.clear();
}

//-----------------------------------------------------------------------------
const char* FEModel::GetUnits() const
{
    if (m_imp->m_units.empty())
        return nullptr;
    else
        return m_imp->m_units.c_str();
}


//! Evaluates all load curves at the specified time
void FEModel::EvaluateLoadControllers(double time)
{
    const int NLC = LoadControllers();
    for (int i = 0; i < NLC; ++i)
        GetLoadController(i)->Evaluate(time);
}

//-----------------------------------------------------------------------------
bool FEModel::EvaluateLoadParameters()
{
    feLog("\n");
    int NLC = LoadControllers();
    for (int i = 0; i < (int)m_imp->m_Param.size(); ++i)
    {
        Impl::LoadParam& pi = m_imp->m_Param[i];
        int nlc = pi.lc;
        if ((nlc >= 0) && (nlc < NLC))
        {
            double s = GetLoadController(nlc)->Value();
            FEParam* p = pi.param;
            FEObjectBase* parent = dynamic_cast<FEObjectBase*>(p->owner());
            if (m_imp->m_printParams)
            {
                if (parent && (parent->GetName().empty() == false))
                {
                    const char* pname = parent->GetName().c_str();
                    feLog("Setting parameter \"%s.%s\" to : ", pname, p->name());
                }
                else
                    feLog("Setting parameter \"%s\" to : ", p->name());
            };
            assert(p->IsVolatile());
            switch (p->type())
            {
                case FE_PARAM_INT:
                {
                    p->value<int>() = (int)s;
                    if (m_imp->m_printParams)
                        feLog("%d\n", p->value<int>());
                }
                break;
                case FE_PARAM_DOUBLE:
                {
                    p->value<double>() = pi.m_scl * s;
                    if (m_imp->m_printParams)
                        feLog("%lg\n", p->value<double>());
                }
                break;
                case FE_PARAM_BOOL:
                {
                    p->value<bool>() = (s > 0 ? true : false);
                    if (m_imp->m_printParams)
                        feLog("%s\n", (p->value<bool>() ? "true" : "false"));
                }
                break;
                case FE_PARAM_VEC3D:
                {
                    Vector3d& v = p->value<Vector3d>();
                    p->value<Vector3d>() = pi.m_vscl * s;
                    if (m_imp->m_printParams)
                        feLog("%lg, %lg, %lg\n", v.x, v.y, v.z);
                }
                break;
                case FE_PARAM_DOUBLE_MAPPED:
                {
                    FEParamDouble& v = p->value<FEParamDouble>();
                    double c = 1.0;
                    if (v.isConst())
                        c = v.constValue();
                    v.SetScaleFactor(s * pi.m_scl);
                    if (m_imp->m_printParams)
                        feLog("%lg\n", c * p->value<FEParamDouble>().GetScaleFactor());
                }
                break;
                case FE_PARAM_VEC3D_MAPPED:
                {
                    FEParamVec3& v = p->value<FEParamVec3>();
                    v.SetScaleFactor(s * pi.m_scl);
                    if (m_imp->m_printParams)
                        feLog("%lg\n", v.GetScaleFactor());
                }
                break;
                default:
                    feLog("\n");
                    assert(false);
            }
        }
        else
        {
            feLogError("Invalid load curve ID");
            return false;
        }
    }
    feLog("\n");

    return true;
}