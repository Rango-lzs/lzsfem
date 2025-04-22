#pragma once
#include "FEStepComponent.h"
#include "FEGlobalVector.h"
#include "femcore/FEDofList.h"

//-----------------------------------------------------------------------------
class FELinearSystem;

//-----------------------------------------------------------------------------
//! This class is the base class for all classes that affect the state of the model
//! and contribute directly to the residual and the global stiffness matrix. This
//! includes most boundary loads, body loads, contact, etc.
class FEM_EXPORT FEModelLoad : public FEStepComponent
{
    DECLARE_META_CLASS(FEModelLoad, FEStepComponent);

public:
	//! constructor
	FEModelLoad(FEModel* pfem);

	const FEDofList& GetDofList() const;
	
	void Serialize(DumpStream& ar) override;

public:
	// all classes derived from this base class must implement
	// the following functions.

	//! evaluate the contribution to the external load vector
	virtual void LoadVector(FEGlobalVector& R);

	//! evaluate the contribution to the global stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS);

protected:
	FEDofList	m_dof;
};
