/*********************************************************************
 * \file   FESolver.h
 * \brief  Refactored FE Solver base class
 *
 * \author Refactored
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "datastructure/Matrix.h"
#include "datastructure/vector_operator.h"
#include "femcore/FEDofList.h"
#include "femcore/FEObjectBase.h"
#include "femcore/FETimeInfo.h"
#include "femcore/Timer.h"
#include "femcore/fecore_enum.h"

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class FEGlobalMatrix;
class LinearSolver;
class FEGlobalVector;

//-----------------------------------------------------------------------------
// Solution variable
class FESolutionVariable
{
public:
    FESolutionVariable(const char* szname, FEDofList* dofs = nullptr, int order = 2)
        : m_szname(szname), m_dofs(dofs), m_order(order)
    {
    }

public:
    FEDofList* m_dofs;     // the dof list
    int m_order;           // the order of interpolation (0 = constant, 1 = linear, 2 = default)
    const char* m_szname;  // name of solution variable
};

//-----------------------------------------------------------------------------
// Structure identifying nodal dof info
struct FEM_EXPORT FENodalDofInfo
{
    int m_eq = -1;              // equation number
    int m_node = -1;            // 0-based index into mesh
    int m_dof = -1;             // index into nodal m_ID array
    const char* szdof = nullptr;
};

//-----------------------------------------------------------------------------
//! Base class for all FE solvers
//! 
//! This class provides the common interface for all solvers.
//! Derived classes implement specific solution strategies.
class FEM_EXPORT FESolver : public FEObjectBase
{
    DECLARE_META_CLASS(FESolver, FEObjectBase);

public:
    //! constructor
    FESolver();

    //! destructor
    virtual ~FESolver();

public:
    //! Data serialization
    void Serialize(DumpStream& ar) override;

    //! Initialize the solver
    virtual bool Init() override;

    //! Clean up solver data
    virtual void Clean();

    //! Rewind the solver (called when time step fails)
    virtual void Rewind();

    //! Reset solver to initial state
    virtual void Reset();

    //! Initialize linear equation system
    virtual bool InitEquations();

    //! Initialize the step
    virtual bool InitStep(double time) = 0;

    //! Solve an analysis step
    virtual bool SolveStep() = 0;

    //! Update the state of the model
    virtual void Update(std::vector<double>& u);

public:
    //! Get the size of a partition
    int GetPartitionSize(int partition);

    //! Get the current stiffness matrix
    virtual FEGlobalMatrix* GetStiffnessMatrix() = 0;

    //! Get the current load vector
    virtual std::vector<double> GetLoadVector() = 0;

    //! Get the linear solver
    virtual LinearSolver* GetLinearSolver() = 0;

    //! Get the solution vector
    virtual std::vector<double> GetSolutionVector() const;

    //! Matrix symmetry flag
    int MatrixSymmetryFlag() const;

    //! Get matrix type
    MatrixType MatType() const;

    //! Build the matrix profile
    virtual void BuildMatrixProfile(FEGlobalMatrix& G, bool breset);

    //! Check if dofs are active
    bool HasActiveDofs(const FEDofList& dof);

    //! Get active dof map
    int GetActiveDofMap(std::vector<int>& dofMap);

    //! Get DOF info from equation number
    FENodalDofInfo GetDOFInfoFromEquation(int ieq);

    //! Extract solution norm for specific dofs
    double ExtractSolutionNorm(const std::vector<double>& v, const FEDofList& dofs) const;

protected:
    //! Add equations
    void AddEquations(int neq, int partition = 0);

    //! Set the linear system partitions
    void SetPartitions(const std::vector<int>& part);

    //! Add a solution variable
    void AddSolutionVariable(FEDofList* dofs, int order, const char* szname);

public:
    // Solver parameters
    bool m_bwopt;               //!< bandwidth optimization flag
    int m_msymm;                //!< matrix symmetry flag
    int m_neq;                  //!< number of equations
    std::vector<int> m_part;    //!< partitions of linear system
    int m_eq_scheme;            //!< equation number scheme (used in InitEquations)
    std::vector<int> m_dofMap;  //!< equation to dof mapping

    // Statistics
    int m_nrhs;     //!< number of RHS evaluations
    int m_niter;    //!< number of iterations
    int m_nref;     //!< number of stiffness reformations (current step)
    int m_ntotref;  //!< total number of stiffness reformations

protected:
    // Solution variables
    std::vector<FESolutionVariable> m_Var;

    DECLARE_PARAM_LIST();
};
