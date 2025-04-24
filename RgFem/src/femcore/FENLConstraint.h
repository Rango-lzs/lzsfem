#pragma once
#include "femcore/Solver/FESolver.h"
#include "FEStepComponent.h"
#include "FEGlobalVector.h"
#include "FEGlobalMatrix.h"
#include "FETimeInfo.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints enforced using an augmented Lagrangian method.

//! The constraint must provide a residual (force) contribution, its stiffness matrix,
//! and an augmentation function.
//!
class FEM_EXPORT FENLConstraint : public FEStepComponent
{
    DECLARE_META_CLASS(FENLConstraint, FEStepComponent);

public:
	FENLConstraint(FEModel* pfem);
	virtual ~FENLConstraint();

	// clone the constraint
	virtual void CopyFrom(FENLConstraint* plc) {}

	// initialize equations
	// Overridden by constraint classes that need to allocate more equations,
	// e.g. for Lagrange Multipliers.
	virtual int InitEquations(int neq) { return 0; }

public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	virtual void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) = 0;

	// Evaluates the contriubtion to the stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) = 0;

	// Performs an augmentation step
	virtual bool Augment(int naug, const FETimeInfo& tp) { return true; }

	// Build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// reset the state data
	virtual void Reset() {}

	// called at start of time step
	virtual void PrepStep() {}

	// update 
	using FEModelComponent::Update;
	virtual void Update(const std::vector<double>& ui) {}
	virtual void Update(const std::vector<double>& Ui, const std::vector<double>& ui) { Update(ui); }
	virtual void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) {}
};
