#pragma once
#include "femcore/fem_export.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class DumpStream;

//-----------------------------------------------------------------------------
// Convenience class for creating a list of degrees of freedom, without the 
// need to go through the FEModel class;

//DofList表示求解的系统开启了那些自由度如：u,v,w,rx,ry,rz
class FEM_EXPORT FEDofList
{
public:
    FEDofList();
	FEDofList(FEModel* fem);
	FEDofList(const FEDofList& dofs);

	// assignment operator
	void operator = (const FEDofList& dofs);
	void operator = (const std::vector<int>& dofs);

	// clear the list
	void Clear();

	// Add a degree of freedom
	bool AddDof(const char* szdof);

	// Add a degree of freedom
	bool AddDof(int ndof);

	// Add all the dofs of a variable
	bool AddVariable(const char* szvar);

	// Add all the dofs a variable
	bool AddVariable(int nvar);

	// Add degrees of freedom
	bool AddDofs(const FEDofList& dofs);

	// is the list empty
	bool IsEmpty() const;

	// number of dofs in the list
	int Size() const;

	// access a dof
	int operator [] (int n) const;

	// serialization
	void Serialize(DumpStream& ar);

	// see if this dof list contains all the dofs of a FEDofList
	bool Contains(int dof);
	bool Contains(const FEDofList& dof);

	// Get the interpolation order 
	int InterpolationOrder(int index) const;

private:
	FEModel*			m_fem;
	std::vector<int>	m_dofList;
};
