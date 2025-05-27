#pragma once
#include "femcore/Solver/FELineSearch.h"
#include "femcore/FETimeInfo.h"
#include "FENewtonStrategy.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
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
struct ConvergenceInfo
{
    int nvar;        // corresponding solution variable
    double tol;      // convergence tolerance
    double norm0;    // initial norm
    double normi;    // current incremental norm
    double norm;     // current total norm
    double maxnorm;  // the max norm

    ConvergenceInfo()
    {
        nvar = -1;
        tol = 0.0;
        norm0 = 0.0;
        normi = 0.0;
        norm = 0.0;
        maxnorm = 0.0;
    }

    bool IsConverged() const
    {
        return (tol > 0 ? normi <= (tol * tol) * norm : true);
    }
};

//-----------------------------------------------------------------------------
// This class defines the base class for Newton-type solvers.
class FEM_EXPORT FENewtonSolver : public FESolver
{
    DECLARE_META_CLASS(FENewtonSolver, FESolver);

public:
    //! constructor
    FENewtonSolver(FEModel* pfem);

    //! destrcutor
    ~FENewtonSolver();

    //! Set the default solution strategy
    void SetDefaultStrategy(QN_STRATEGY qn);

    //! Check the zero diagonal
    void CheckZeroDiagonal(bool bcheck, double ztol = 0.0);

public:  // overloaded from FESolver
    //! Initialization
    bool Init() override;

    //! return number of equations
    int NumberOfEquations() const
    {
        return m_neq;
    }

    //! Clean up
    void Clean() override;

    //! serialization
    void Serialize(DumpStream& ar) override;

    //! Solve an analysis step
    bool SolveStep() override;

    //! rewind solver
    void Rewind() override;

    //! prep the solver for the QN updates
    virtual void PrepStep();

public:  // Quasi-Newton methods
    //! call this at the start of the quasi-newton loop
    bool QNInit();

    //! Do a qn update
    bool QNUpdate();

    //! solve the equations using QN method (returns line search size)
    double QNSolve();

    //! Force a stiffness reformation during next update
    void QNForceReform(bool b);

    // return line search
    FELineSearch* GetLineSearch();

public:
    //! return the stiffness Matrix
    FEGlobalMatrix* GetStiffnessMatrix() override;

    //! reform the stiffness Matrix
    bool ReformStiffness();

    //! recalculates the shape of the stiffness Matrix
    bool CreateStiffness(bool breset);

    //! get the RHS
    std::vector<double> GetLoadVector() override;

    //! Get the total solution vector (for current Newton iteration)
    virtual void GetSolutionVector(std::vector<double>& U);

    //! Get the total solution vector
    std::vector<double> GetSolutionVector() const override;

public:
    //! do augmentations
    bool DoAugmentations();

    //! solve the equations
    void SolveEquations(std::vector<double>& u, std::vector<double>& R);

    //! do a line search
    double DoLineSearch();

public:
    //! Set the solution strategy
    void SetSolutionStrategy(FENewtonStrategy* pstrategy);

    //! solve the linear system of equations
    void SolveLinearSystem(std::vector<double>& x, std::vector<double>& R);

    //! Do a Quasi-Newton step
    //! This is called from SolveStep.
    virtual bool Quasin();

    //! calculates the global stiffness Matrix (needs to be overwritten by derived classes)
    virtual bool StiffnessMatrix();

    //! this is the new method for building the stiffness Matrix
    virtual bool StiffnessMatrix(FELinearSystem& LS);

    //! calculates the global residual vector (needs to be overwritten by derived classes)
    virtual bool Residual(std::vector<double>& R) = 0;

    //! Check convergence. Derived classes that don't override Quasin, should implement this
    //! niter = iteration number
    //! ui    = search direction
    //! ls    = line search factor
    virtual bool CheckConvergence(int niter, const std::vector<double>& ui, double ls);

    //! return the linear solver
    LinearSolver* GetLinearSolver() override;

    //! Add a solution variable from a doflist
    void AddSolutionVariable(FEDofList* dofs, int order, const char* szname, double tol);

    //! Update the state of the model
    void Update(std::vector<double>& u) override;

    //! TODO: This is a helper function to get around an issue with the current implementation
    //        regarding prescribed displacements. The purpose of this Update2 function is to update
    //        all degrees of freedom, including prescribed ones. This is currently only used by the JFNKMatrix class.
    //        and overridden in FESolidSolver2.
    virtual void Update2(const std::vector<double>& ui);

    //! Update the model
    virtual void UpdateModel();

protected:
    bool AllocateLinearSystem();

public:
    // line search options
    FELineSearch* m_lineSearch;

    // solver parameters
    int m_maxref;           //!< max nr of reformations per time step
    int m_force_partition;  //!< Force a partition of the global Matrix (e.g. for testing with BIPN solver)
    double m_Rtol;          //!< residual convergence norm
    double m_Etol;          //!< energy convergence norm
    double m_Rmin;          //!< min residual value
    double m_Rmax;          //!< max residual value

    // solution strategy
    FENewtonStrategy* m_qnstrategy;  //!< class handling the specific stiffness update logic
    bool m_breformtimestep;          //!< reform at start of time step
    bool m_breformAugment;           //!< reform after each (failed) augmentations
    bool m_bforceReform;             //!< forces a reform in QNInit
    bool m_bdivreform;               //!< reform when diverging
    bool m_bdoreforms;               //!< do reformations

    // counters
    int m_nref;  //!< nr of stiffness retormations

    // Error handling
    bool m_bzero_diagonal;  //!< check for zero diagonals
    double m_zero_tol;      //!< tolerance for zero diagonal

    // linear solver data
    LinearSolver* m_plinsolve;  //!< the linear solver
    FEGlobalMatrix* m_pK;       //!< global stiffness Matrix
    bool m_breshape;            //!< Matrix reshape flag
    bool
        m_persistMatrix;  //!< Don't delete stiffness Matrix until necessary (if true, K is deleted at end of time step)

    // data used by Quasin
    std::vector<double> m_R0;  //!< residual at iteration i-1
    std::vector<double> m_R1;  //!< residual at iteration i
    std::vector<double> m_ui;  //!< solution increment vector
    std::vector<double> m_Ut;  //!< total solution vector
    std::vector<double> m_Ui;  //!< total solution increments of current time step
    std::vector<double> m_up;  //!< solution increment of previous iteration
    std::vector<double> m_Fd;  //!< residual correction due to prescribed degrees of freedom

private:
    double m_ls;  //!< line search factor calculated in last call to QNSolve

private:
    ConvergenceInfo m_residuNorm;                 // residual convergence info
    ConvergenceInfo m_energyNorm;                 // energy convergence info
    std::vector<ConvergenceInfo> m_solutionNorm;  // converge info for solution variables

    DECLARE_PARAM_LIST();
};
