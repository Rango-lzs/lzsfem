#pragma once
#include "femcore/FELinearConstraint.h"
#include "datastructure/table.h"

class FEGlobalMatrix;
class Matrix;

//-----------------------------------------------------------------------------
// This class helps manage all the linear constraints
class FEM_EXPORT FELinearConstraintManager
{
public:
	FELinearConstraintManager(FEModel* fem);
	~FELinearConstraintManager();

	// Clear all constraints
	void Clear();

	// copy data
	void CopyFrom(const FELinearConstraintManager& lcm);

	// serialize linear constraints
	void Serialize(DumpStream& ar);

	// build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& G);

	// add the linear constraint
	void AddLinearConstraint(FELinearConstraint* lc);

	// return number of linear constraints
	int LinearConstraints() const;

	// return linear constraint
	const FELinearConstraint& LinearConstraint(int i) const;

	// return linear constraint
	FELinearConstraint& LinearConstraint(int i);

	//! remove a linear constraint
	void RemoveLinearConstraint(int i);

public:
	// one-time initialization
	bool Initialize();

	// activation
	bool Activate();

	// assemble element residual into global residual
	void AssembleResidual(std::vector<double>& R, std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe);

	// assemble element matrix into (reduced) global matrix
	void AssembleStiffness(FEGlobalMatrix& K, std::vector<double>& R, std::vector<double>& ui, const std::vector<int>& en, const std::vector<int>& lmi, const std::vector<int>& lmj, const Matrix& ke);

	// called before the first reformation for each time step
	void PrepStep();

	// update nodal variables
	void Update();

protected:
	void InitTable();

private:
	FEModel* m_fem;
	std::vector<FELinearConstraint*>	m_LinC;		//!< linear constraints data
	Table2d<int>					m_LCT;		//!< linear constraint table
	std::vector<double>				m_up;		//!< the inhomogenous component of the linear constraint
};
