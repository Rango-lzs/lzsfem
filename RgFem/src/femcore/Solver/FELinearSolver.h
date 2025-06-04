#pragma once
#include "FESolver.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEGlobalVector;
class FEGlobalMatrix;
class FELinearSystem;
class LinearSolver;

//-----------------------------------------------------------------------------
//! Abstract Base class for finite element solution algorithms (i.e. "FE solvers") that require the solution
//! of a linear system of equations.
class FEM_EXPORT FELinearSolver : public FESolver
{
public:
	//! constructor
	FELinearSolver();

	//! Set the degrees of freedom
	void SetDOF(std::vector<int>& dof);

	//! Get the number of equations
	int NumberOfEquations() const;

	//! Get the linear solver
	LinearSolver* GetLinearSolver() override;

public: // from FESolver

	//! solve the step
	bool SolveStep() override;

	//! Initialize and allocate data
	bool Init() override;

	//! Clean up data
	void Clean() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

public: // these functions need to be implemented by the derived class

	//! Evaluate the right-hand side "force" vector
	virtual void ForceVector(FEGlobalVector& R);

	//! Evaluate the stiffness Matrix
	virtual bool StiffnessMatrix(FELinearSystem& K);

	//! Update the model state
	virtual void Update(std::vector<double>& u) override;

protected: // some helper functions

	//! Reform the stiffness Matrix
	bool ReformStiffness();

	//! Create and evaluate the stiffness Matrix
	bool CreateStiffness();

	//! get the stiffness Matrix
	FEGlobalMatrix* GetStiffnessMatrix() override;

	//! get the RHS
	std::vector<double>	GetLoadVector() override;

protected:
	std::vector<double>		m_R;	//!< RHS vector
	std::vector<double>		m_u;	//!< vector containing prescribed values

private:
	LinearSolver*		m_pls;		//!< The linear equation solver
	FEGlobalMatrix*		m_pK;		//!< The global stiffness Matrix

	std::vector<int>		m_dof;	//!< list of active degrees of freedom
	bool			m_breform;	//!< Matrix reformation flag
};
