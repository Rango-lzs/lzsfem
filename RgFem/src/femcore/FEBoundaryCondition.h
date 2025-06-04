#pragma once
#include "femcore/FEStepComponent.h"
#include "FENodeSet.h"
#include "FEDofList.h"

//-----------------------------------------------------------------------------
class FEFacetSet;

//-----------------------------------------------------------------------------
//! This class is the base class of boundary conditions.

//! Boundary conditions set the "bc" state of nodes. The bc-state determines
//! whether or not the dofs of the node will be assigned an equation number. 
//! Currently, there are two boundary conditions: a fixed (FEFixedBC) and a
//! prescribed (FEPrescribedBC) boundary condition. 
class FEM_EXPORT FEBoundaryCondition : public FEStepComponent
{
    DECLARE_META_CLASS(FEBoundaryCondition, FEStepComponent);

public:
	//! constructor
	FEBoundaryCondition();

	//! desctructor
	~FEBoundaryCondition();

	//! fill the prescribed values
	virtual void PrepStep(std::vector<double>& u, bool brel = true);

	// copy data from another class
	virtual void CopyFrom(FEBoundaryCondition* pbc) = 0;
    
    // repair BC if needed
    virtual void Repair() {}

	void Serialize(DumpStream& ar) override;

	// TODO: Temporary construction to update some special boundary conditions in FEModel::Update
	//       Will likely remove this at some point.
	virtual void UpdateModel() {}


public:
	// set the dof list
	void SetDOFList(int ndof);
	void SetDOFList(const std::vector<int>& dofs);
	void SetDOFList(const FEDofList& dofs);

	const FEDofList& GetDofList() const { return m_dof; }

protected:
	FEDofList	m_dof;	// the dof list for the BC
};
