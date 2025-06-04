#pragma once
#include "FEStepComponent.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
// TODO: I like to inherit this from FENLConstraint and potentially eliminate
//       The distinction between a nonlinear constraint and a contact interface.
//       Since a contact interface essentially is a nonlinear constraint, I think
//       this may make things a lot easier. I already made the function definitions consistent
//       but am hesitant to push this through at this point. 
class FEM_EXPORT FESurfacePairConstraint : public FEStepComponent
{
    DECLARE_META_CLASS(FESurfacePairConstraint, FEStepComponent);

public:
	//! constructor
	FESurfacePairConstraint();

public:
	//! return the primary surface
	virtual FESurface* GetPrimarySurface() = 0;

	//! return the secondary surface
	virtual FESurface* GetSecondarySurface () = 0;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	virtual bool UseNodalIntegration() = 0;

	//! create a copy of this interface
	virtual void CopyFrom(FESurfacePairConstraint* pci) {}

public:
	// Build the Matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// reset the state data
	virtual void Reset() {}

	// do the augmentation
	virtual bool Augment(int naug, const FETimeInfo& tp) { return true; }

	// allocate equations for lagrange multipliers
	// (should return the number of equations to be allocated)
	virtual int InitEquations(int neq);

	// update based on solution (use for updating Lagrange Multipliers)
	virtual void Update(std::vector<double>& ui);

	using FEModelComponent::Update;
};
