/*********************************************************************
 * \file   FENewtonSolver.h
 * \brief  Refactored Newton solver base class
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "femcore/NewtonSolver/SolverBase.h"
#include "femcore/Solver/FELineSearch.h"
#include "femcore/Solver/FENewtonStrategy.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEGlobalMatrix;
class FELinearSystem;

//-----------------------------------------------------------------------------
enum QN_STRATEGY
{
    QN_BFGS,
    QN_BROYDEN,
    QN_JFNK
};

//-----------------------------------------------------------------------------
//! Convergence information structure
struct ConvergenceInfo
{
    int nvar;        // corresponding solution variable
    double tol;      // convergence tolerance
    double norm0;    // initial norm
    double normi;    // current incremental norm
    double norm;     // current total norm
    double maxnorm;  // maximum norm

    ConvergenceInfo()
        : nvar(-1), tol(0.0), norm0(0.0), normi(0.0), norm(0.0), maxnorm(0.0)
    {
    }

    bool IsConverged() const
    {
        return (tol > 0 ? normi <= (tol * tol) * norm : true);
    }
};

//-----------------------------------------------------------------------------
//! Base class for Newton-type solvers
//! 
//! This class implements the common Newton-Raphson iteration logic.
//! Derived classes implement specific physics (static, implicit dynamic, etc.)
class FEM_EXPORT FENewtonSolver : public FESolver
{
    DECLARE_META_CLASS(FENewtonSolver, FESolver);

public:
    //! Constructor
    FENewtonSolver();

    //! Destructor
    virtual ~FENewtonSolver();

public:
    //! Initialization
    bool Init() override;

    //! Clean up
    void Clean() override;

    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! Solve an analysis step
    bool SolveStep() override;

    //! Rewind solver
    void Rewind() override;

public:  // Newton iteration methods
    //! Prepare for the QN updates
    virtual void PrepStep();

    //! Initialize quasi-Newton iteration
    bool QNInit();

    //! Perform a quasi-Newton update
    bool QNUpdate();

    //! Solve equations using QN method (returns line search size)
    double QNSolve();

    //! Perform the Quasi-Newton iteration loop
    virtual bool Quasin();

    //! Force stiffness reformation during next update
    void QNForceReform(bool b);

public:  // Stiffness matrix methods
    //! Get the stiffness matrix
    FEGlobalMatrix* GetStiffnessMatrix() override;

    //! Reform the stiffness matrix
    bool ReformStiffness();

    //! Recalculate the matrix profile
    bool CreateStiffness(bool breset);

    //! Calculate the global stiffness matrix
    virtual bool StiffnessMatrix();

    //! Build stiffness matrix (new method using linear system)
    virtual bool StiffnessMatrix(FELinearSystem& LS);

public:  // Residual methods
    //! Calculate the global residual vector
    virtual bool Residual(std::vector<double>& R) = 0;

    //! Get the RHS vector
    std::vector<double> GetLoadVector() override;

    //! Get the solution vector
    std::vector<double> GetSolutionVector() const override;

    //! Get the total solution vector for current iteration
    virtual void GetSolutionVector(std::vector<double>& U);

public:  // Update methods
    //! Update the model state
    void Update(std::vector<double>& u) override;

    //! Update the model (for use with prescribed DOFs)
    virtual void Update2(const std::vector<double>& ui);

    //! Update the model after convergence
    virtual void UpdateModel();

public:  // Convergence and augmentation
    //! Check convergence
    virtual bool CheckConvergence(int niter, const std::vector<double>& ui, double ls);

    //! Perform augmentations
    bool DoAugmentations();

public:  // Linear solver methods
    //! Get the linear solver
    LinearSolver* GetLinearSolver() override;

    //! Solve the linear system
    void SolveLinearSystem(std::vector<double>& x, std::vector<double>& R);

    //! Solve the equations
    void SolveEquations(std::vector<double>& u, std::vector<double>& R);

public:  // Line search
    //! Perform line search
    double DoLineSearch();

    //! Get line search object
    FELineSearch* GetLineSearch();

public:  // Configuration methods
    //! Set the default solution strategy
    void SetDefaultStrategy(QN_STRATEGY qn);

    //! Set the solution strategy
    void SetSolutionStrategy(FENewtonStrategy* pstrategy);

    //! Enable zero diagonal checking
    void CheckZeroDiagonal(bool bcheck, double ztol = 0.0);

    //! Add a solution variable
    void AddSolutionVariable(FEDofList* dofs, int order, const char* szname, double tol);

    //! Get equation count
    int NumberOfEquations() const { return m_neq; }

protected:
    //! Allocate linear system
    bool AllocateLinearSystem();

public:
    // Solver parameters
    int m_maxref;           //!< max reformations per time step
    int m_force_partition;  //!< force matrix partition
    double m_Rtol;          //!< residual tolerance
    double m_Etol;          //!< energy tolerance
    double m_Rmin;          //!< minimum residual value
    double m_Rmax;          //!< maximum residual value

    // Solution strategy
    FENewtonStrategy* m_qnstrategy;  //!< QN strategy
    bool m_breformtimestep;          //!< reform at start of time step
    bool m_breformAugment;           //!< reform after failed augmentations
    bool m_bforceReform;             //!< force reform in QNInit
    bool m_bdivreform;               //!< reform when diverging
    bool m_bdoreforms;               //!< enable reformations

    // Line search
    FELineSearch* m_lineSearch;

    // Error handling
    bool m_bzero_diagonal;  //!< check for zero diagonals
    double m_zero_tol;      //!< tolerance for zero diagonal

//protected:
    // Linear solver data
    LinearSolver* m_plinsolve;  //!< linear solver
    FEGlobalMatrix* m_pK;       //!< global stiffness matrix
    bool m_breshape;            //!< matrix reshape flag
    bool m_persistMatrix;       //!< persist matrix until necessary

    // Iteration data
    std::vector<double> m_R0;  //!< residual at iteration i-1
    std::vector<double> m_R1;  //!< residual at iteration i
    std::vector<double> m_ui;  //!< solution increment vector
    std::vector<double> m_Ut;  //!< total solution vector
    std::vector<double> m_Ui;  //!< total solution increments
    std::vector<double> m_up;  //!< previous iteration increment
    std::vector<double> m_Fd;  //!< prescribed DOF correction

    double m_ls;  //!< line search factor

    // Convergence information
    ConvergenceInfo m_residuNorm;                 //!< residual convergence
    ConvergenceInfo m_energyNorm;                 //!< energy convergence
    std::vector<ConvergenceInfo> m_solutionNorm;  //!< solution convergence

    DECLARE_PARAM_LIST();
};
