#pragma once
#include "FENodalBC.h"

//-----------------------------------------------------------------------------
//! This class represents a fixed degree of freedom
//! This boundary conditions sets the BC attribute of the nodes in the nodeset
//! to DOF_FIXED when activated.
class FEM_EXPORT FEFixedBC : public FENodalBC
{
public:
	//! constructors
	FEFixedBC(FEModel* pfem);
	FEFixedBC(FEModel* pfem, int dof, FENodeSet* ps);

	//! initialization
	bool Init() override;

	//! activation
	void Activate() override;

	//! deactivation
	void Deactivate() override;

	void CopyFrom(FEBoundaryCondition* bc) override;
};

//-----------------------------------------------------------------------------
// This class is obsolete, but provides a direct parameterization of the base class.
// This is maintained for backward compatibility with older feb files.
class FEM_EXPORT FEFixedDOF : public FEFixedBC
{
	DECLARE_META_CLASS(FEFixedDOF, FEFixedBC);

public:
	FEFixedDOF(FEModel* fem);

	void SetDOFS(const std::vector<int>& dofs);

	bool Init() override;

protected:
	std::vector<int>	m_dofs;
	DECLARE_PARAM_LIST();
};
