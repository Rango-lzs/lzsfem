#pragma once
#include "femcore/FEMeshPartition.h"

// forward declaration of material class
class FEMaterial;

// Base class for solid and shell parts. Domains can also have materials assigned.
class FEM_EXPORT FEDomain : public FEMeshPartition
{
    DECLARE_META_CLASS(FEDomain, FEMeshPartition);

public:
	FEDomain(int nclass, FEModel* fem);

	//! get the material of this domain
	virtual FEMaterial* GetMaterial() { return 0; }

	// assign a material to this domain
	virtual void SetMaterial(FEMaterial* pm);

	//! set the material ID of all elements
	void SetMatID(int mid);

	//! Allocate material point data for the elements
	//! This is called after elements get read in from the input file.
	//! And must be called before material point data can be accessed.
	//! \todo Perhaps I can make this part of the "creation" routine
	void CreateMaterialPointData();

	// serialization
	void Serialize(DumpStream& ar) override;

	//! augmentation
	// NOTE: This is here so that the FESolver can do the augmentations
	// for the 3-field hex/shell domains.
	virtual bool Augment(int naug) { return true; }

	// create function
	virtual bool Create(int elements, FE_Element_Spec espec) = 0;

public:
	//! Get the list of dofs on this domain
	virtual const FEDofList& GetDOFList() const = 0;

	//! Unpack the LM data for an element of this domain
	virtual void UnpackLM(FEElement& el, std::vector<int>& lm);

	//! build the Matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M);

	//! Activate the domain
	virtual void Activate();

protected:
	// helper function for activating dof lists
	void Activate(const FEDofList& dof);

	// helper function for unpacking element dofs
	void UnpackLM(FEElement& el, const FEDofList& dof, std::vector<int>& lm);
};
