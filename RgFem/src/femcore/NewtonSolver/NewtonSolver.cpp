/*********************************************************************
 * \file   FENewtonSolver.cpp
 * \brief  Implementation of Newton solver base class
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#include "femcore/NewtonSolver/NewtonSolver.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/Matrix/FEGlobalMatrix.h"
#include "femcore/FELinearSystem.h"
#include "femcore/Solver/LinearSolver.h"
#include "logger/log.h"

DEFINE_META_CLASS(FENewtonSolver, FESolver,"newton");

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FENewtonSolver, FESolver)
    ADD_PARAMETER(m_maxref, "max_refs");
    ADD_PARAMETER(m_Rtol, "rtol");
    ADD_PARAMETER(m_Etol, "etol");
    ADD_PARAMETER(m_Rmin, "rmin");
    ADD_PARAMETER(m_Rmax, "rmax");
    ADD_PARAMETER(m_breformtimestep, "reform_each_time_step");
    ADD_PARAMETER(m_bdivreform, "diverge_reform");
    ADD_PARAMETER(m_breformAugment, "reform_augment");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FENewtonSolver::FENewtonSolver()
{
    // Solver parameters
    m_maxref = 15;
    m_force_partition = 0;
    m_Rtol = 1e-5;
    m_Etol = 1e-4;
    m_Rmin = 1e-20;
    m_Rmax = 1e20;
    
    // Strategy
    m_qnstrategy = nullptr;
    m_breformtimestep = true;
    m_breformAugment = false;
    m_bforceReform = false;
    m_bdivreform = false;
    m_bdoreforms = true;
    
    // Line search
    m_lineSearch = new FELineSearch(this);
    
    // Error handling
    m_bzero_diagonal = false;
    m_zero_tol = 0.0;
    
    // Linear solver
    m_plinsolve = nullptr;
    m_pK = nullptr;
    m_breshape = true;
    m_persistMatrix = false;
    
    // Counters
    m_nref = 0;
    
    m_ls = 1.0;
}

//-----------------------------------------------------------------------------
FENewtonSolver::~FENewtonSolver()
{
    if (m_lineSearch) delete m_lineSearch;
    if (m_qnstrategy) delete m_qnstrategy;
    if (m_pK) delete m_pK;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::Init()
{
    if (!FESolver::Init()) return false;
    
    // Allocate linear system
    if (!AllocateLinearSystem()) return false;
    
    // Initialize line search
    //if (m_lineSearch && !m_lineSearch->Init(this))
    //    return false;
    
    // Initialize QN strategy
    if (m_qnstrategy && !m_qnstrategy->Init())
        return false;
    
    return true;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Clean()
{
    FESolver::Clean();
    
    m_R0.clear();
    m_R1.clear();
    m_ui.clear();
    m_Ut.clear();
    m_Ui.clear();
    m_up.clear();
    m_Fd.clear();
    
    if (m_pK && !m_persistMatrix)
    {
        delete m_pK;
        m_pK = nullptr;
    }
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Serialize(DumpStream& ar)
{
    FESolver::Serialize(ar);
    
    ar & m_nref;
    ar & m_R0 & m_R1;
    ar & m_ui & m_Ut & m_Ui & m_up & m_Fd;
    
    if (m_qnstrategy)
        m_qnstrategy->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Rewind()
{
    FESolver::Rewind();
    
    // Force reformation on next step
    m_bforceReform = true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::AllocateLinearSystem()
{
    FEModel* fem = GetFEModel();
    if (!fem) return false;
    
    // Create stiffness matrix if needed
    if (!m_pK)
    {
        // Create the appropriate matrix type
        //FECoreKernel& fecore = FECoreKernel::GetInstance();
        //m_pK = fecore.CreateGlobalMatrix(MatType());
        
        if (!m_pK)
        {
            RgLogError("Failed to create stiffness matrix");
            return false;
        }
    }
    
    // Allocate vectors
    m_R0.assign(m_neq, 0.0);
    m_R1.assign(m_neq, 0.0);
    m_ui.assign(m_neq, 0.0);
    m_Ut.assign(m_neq, 0.0);
    m_Ui.assign(m_neq, 0.0);
    m_up.assign(m_neq, 0.0);
    m_Fd.assign(m_neq, 0.0);
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::SolveStep()
{
    FEModel* fem = GetFEModel();
    
    // Initialize counters
    m_niter = 0;
    m_nref = 0;
    m_nrhs = 0;
    
    // Prepare step
    PrepStep();
    
    // Perform Newton iteration
    bool converged = Quasin();
    
    if (converged)
    {
        // Update model
        UpdateModel();
    }
    
    return converged;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::PrepStep()
{
    // Zero total displacement increment
    std::fill(m_Ui.begin(), m_Ui.end(), 0.0);
    
    // Reset reform flag
    if (m_breformtimestep)
        m_bforceReform = true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::Quasin()
{
    FEModel* fem = GetFEModel();
    
    // Initialize QN iteration
    if (!QNInit()) return false;
    
    // Perform QN iterations
    bool converged = false;
    do
    {
        // Perform QN update
        if (!QNUpdate()) return false;
        
        // Solve equations
        double s = QNSolve();
        
        // Check convergence
        converged = CheckConvergence(m_niter, m_ui, s);
        
        // Increment iteration counter
        m_niter++;
        
        // Check max iterations
        if (m_niter >= m_maxref * 10)  // Some reasonable limit
        {
            RgLogError("Max iterations exceeded");
            return false;
        }
        
    } while (!converged);
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::QNInit()
{
    // Calculate residual
    m_nrhs++;
    if (!Residual(m_R0)) return false;
    
    // Reform stiffness if needed
    bool reform = m_bforceReform || m_breformtimestep;
    
    if (reform)
    {
        if (!ReformStiffness()) return false;
        m_bforceReform = false;
    }
    
    // Initialize strategy
    //if (m_qnstrategy && !m_qnstrategy->InitStep())
    //    return false;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::QNUpdate()
{
    // Calculate residual
    m_nrhs++;
    if (!Residual(m_R1)) return false;
    
    // Check if we need to reform
    bool reform = false;
    
    if (m_qnstrategy)
    {
        //reform = m_qnstrategy->RequireReform(m_niter, m_ui, m_R1);
    }
    
    if (reform && m_bdoreforms)
    {
        if (!ReformStiffness()) return false;
    }
    
    // Update strategy
    if (m_qnstrategy)
    {
        //if (!m_qnstrategy->Update(m_ui, m_R0, m_R1))
       //     return false;
    }
    
    // Copy residual
    m_R0 = m_R1;
    
    return true;
}

//-----------------------------------------------------------------------------
double FENewtonSolver::QNSolve()
{
    // Solve linear system
    SolveLinearSystem(m_ui, m_R1);
    
    // Perform line search
    double s = 1.0;
    if (m_lineSearch && m_lineSearch->m_LSiter > 0)
    {
        s = DoLineSearch();
    }
    
    // Update solution
    for (int i = 0; i < m_neq; ++i)
    {
        m_Ui[i] += s * m_ui[i];
        m_Ut[i] += s * m_ui[i];
    }
    
    // Update model
    Update(m_ui);
    
    m_ls = s;
    return s;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::QNForceReform(bool b)
{
    m_bforceReform = b;
}

//-----------------------------------------------------------------------------
FELineSearch* FENewtonSolver::GetLineSearch()
{
    return m_lineSearch;
}

//-----------------------------------------------------------------------------
FEGlobalMatrix* FENewtonSolver::GetStiffnessMatrix()
{
    return m_pK;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::ReformStiffness()
{
    // Build stiffness matrix
    if (!StiffnessMatrix()) return false;
    
    m_nref++;
    m_ntotref++;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::CreateStiffness(bool breset)
{
    if (!m_pK) return false;
    
    // Build matrix profile
    BuildMatrixProfile(*m_pK, breset);
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::StiffnessMatrix()
{
    // Create linear system
    FELinearSystem LS(this, *m_pK, m_Fd, m_ui, m_breshape);
    
    // Call new interface
    if (!StiffnessMatrix(LS)) return false;
    
    m_breshape = false;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::StiffnessMatrix(FELinearSystem& LS)
{
    // Default implementation - derived classes should override
    return false;
}

//-----------------------------------------------------------------------------
std::vector<double> FENewtonSolver::GetLoadVector()
{
    return m_R1;
}

//-----------------------------------------------------------------------------
std::vector<double> FENewtonSolver::GetSolutionVector() const
{
    return m_Ut;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::GetSolutionVector(std::vector<double>& U)
{
    U = m_Ut;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::DoAugmentations()
{
    FEModel* fem = GetFEModel();
    
    // Perform augmentations
    bool converged = false; // Augment();
    
    if (!converged && m_breformAugment)
    {
        m_bforceReform = true;
    }
    
    return converged;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::SolveEquations(std::vector<double>& u, std::vector<double>& R)
{
    SolveLinearSystem(u, R);
}

//-----------------------------------------------------------------------------
double FENewtonSolver::DoLineSearch()
{
    if (!m_lineSearch) return 1.0;
    
    //return m_lineSearch->DoLineSearch(m_ui, m_R1);
}

//-----------------------------------------------------------------------------
void FENewtonSolver::SetSolutionStrategy(FENewtonStrategy* pstrategy)
{
    if (m_qnstrategy) delete m_qnstrategy;
    m_qnstrategy = pstrategy;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::SolveLinearSystem(std::vector<double>& x, std::vector<double>& R)
{
    if (!m_plinsolve || !m_pK)
    {
        RgLogError("Linear solver not initialized");
        return;
    }
    
    // Solve K * x = R
    if (!m_plinsolve->BackSolve(x, R))
    {
        RgLogError("Linear solver failed");
    }
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::CheckConvergence(int niter, const std::vector<double>& ui, double ls)
{
    // Calculate norms
    m_residuNorm.normi = 0.0;
    m_energyNorm.normi = 0.0;
    
    for (int i = 0; i < m_neq; ++i)
    {
        m_residuNorm.normi += m_R1[i] * m_R1[i];
        m_energyNorm.normi += m_R1[i] * ui[i];
    }
    
    if (niter == 0)
    {
        m_residuNorm.norm0 = m_residuNorm.normi;
        m_energyNorm.norm0 = m_energyNorm.normi;
    }
    
    // Set tolerances
    m_residuNorm.tol = m_Rtol;
    m_energyNorm.tol = m_Etol;
    
    // Check convergence
    bool rconv = m_residuNorm.IsConverged();
    bool econv = m_energyNorm.IsConverged();
    
    return rconv && econv;
}

//-----------------------------------------------------------------------------
LinearSolver* FENewtonSolver::GetLinearSolver()
{
    return m_plinsolve;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::AddSolutionVariable(FEDofList* dofs, int order, const char* szname, double tol)
{
    FESolver::AddSolutionVariable(dofs, order, szname);
    
    // Add convergence info
    ConvergenceInfo ci;
    ci.nvar = (int)m_Var.size() - 1;
    ci.tol = tol;
    m_solutionNorm.push_back(ci);
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Update(std::vector<double>& u)
{
    // Default implementation
    FESolver::Update(u);
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Update2(const std::vector<double>& ui)
{
    // Default - same as Update
    std::vector<double> u = ui;
    Update(u);
}

//-----------------------------------------------------------------------------
void FENewtonSolver::UpdateModel()
{
    FEModel* fem = GetFEModel();
    if (fem) fem->Update();
}

//-----------------------------------------------------------------------------
void FENewtonSolver::SetDefaultStrategy(QN_STRATEGY qn)
{
    if (m_qnstrategy) delete m_qnstrategy;
    
    switch (qn)
    {
        /*case QN_BFGS:
            m_qnstrategy = new FEBFGSStrategy(this);
            break;
        case QN_BROYDEN:
            m_qnstrategy = new FEBroydenStrategy(this);
            break;
        case QN_JFNK:
            m_qnstrategy = new FEJFNKStrategy(this);
            break;*/
        default:
            m_qnstrategy = nullptr;
    }
}

//-----------------------------------------------------------------------------
void FENewtonSolver::CheckZeroDiagonal(bool bcheck, double ztol)
{
    m_bzero_diagonal = bcheck;
    m_zero_tol = ztol;
}
