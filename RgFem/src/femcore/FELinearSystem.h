#pragma once
#include "femcore/Matrix/FEGlobalMatrix.h"
#include "datastructure/Matrix.h"
#include "femcore/fem_export.h"
#include <vector>

class FESolver;

//-----------------------------------------------------------------------------
// Experimental class to see if all the assembly operations can be moved to a class
// and out of the solver class.
class FEM_EXPORT FELinearSystem
{
public:
	// Constructor
	// Takes a FEGlobalMatrix class K that will store the actual stiffness Matrix
	// and a vector F which contains the assembled contribution of the prescribed 
	// degrees of freedom. The F vector must be added to the "force" vector. The u 
	// vector contains the nodal values of the prescribed degrees of freedom.
	FELinearSystem(FESolver* solver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm);

	// virtual destructor
	virtual ~FELinearSystem();

	// get symmetry flag
	bool IsSymmetric() const;

	// Get the solver that is using this linear system
	FESolver* GetSolver();

public:
	// Assembly routine
	// This assembles the element stiffness Matrix ke into the global Matrix.
	// The contributions of prescribed degrees of freedom will be stored in m_F
	virtual void Assemble(const FEElementMatrix& ke);

	// This assembles a Matrix to the RHS by pre-multiplying the Matrix with the 
	// prescribed value array U and then adding it to F
	void AssembleRHS(std::vector<int>& lm, Matrix& ke, std::vector<double>& U);

	// This assembles a vetor to the RHS
	void AssembleRHS(std::vector<int>& lm, std::vector<double>& fe);

protected:
	bool					m_bsymm;	//!< symmetry flag
	FESolver*				m_solver;
	FEGlobalMatrix&			m_K;	//!< The global stiffness Matrix
	std::vector<double>&	m_F;	//!< Contributions from prescribed degrees of freedom
	std::vector<double>&	m_u;	//!< the array with prescribed values
};
