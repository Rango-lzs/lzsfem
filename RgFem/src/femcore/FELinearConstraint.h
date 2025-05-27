#pragma once
#include "femcore/FEBoundaryCondition.h"
#include <vector>

//-----------------------------------------------------------------------------
//! linear constraint
class FEM_EXPORT FELinearConstraintDOF : public FEObjectBase
{
	DECLARE_META_CLASS(FELinearConstraintDOF, FEObjectBase);
public:
	FELinearConstraintDOF(FEModel* fem);

	public:
		int		node;	// node number
		int		dof;	// degree of freedom
		double	val;	// coefficient value (ignored for parent dof)

private:
	FELinearConstraintDOF(const FELinearConstraintDOF&);
	void operator = (const FELinearConstraintDOF&) {}

	DECLARE_PARAM_LIST();
};

class FEM_EXPORT FELinearConstraint : public FEBoundaryCondition
{
    DECLARE_META_CLASS(FELinearConstraint, FEBoundaryCondition);

public:
	typedef std::vector<FELinearConstraintDOF*>::iterator dof_iterator;

public:
	// constructors
	FELinearConstraint();
	FELinearConstraint(FEModel* pfem);
	FELinearConstraint(const FELinearConstraint& LC);

	~FELinearConstraint();

	void Clear();

	// copy data
	void CopyFrom(FEBoundaryCondition* pbc) override;

	// serialization
	void Serialize(DumpStream& ar);

	// initialize the linear constraint
	bool Init();

	// make the constraint active
	void Activate();
	void Deactivate();

	// set the parent degree of freedom
	void SetParentDof(int dof, int node);
	void SetParentNode(int node);
	void SetParentDof(int dof);

	// get the parent dof
	int GetParentDof() const;
	int GetParentNode() const;

	// add a child degree of freedom
	void AddChildDof(int dof, int node, double v);
	void AddChildDof(FELinearConstraintDOF* dof);

	// set the linear constraint offset
	void SetOffset(double d) { m_off = d; }

	// return offset
	double GetOffset() const;

	// get the child DOF
	const FELinearConstraintDOF& GetChildDof(int n) const;

	size_t Size() const;

	dof_iterator begin();

protected:
	FELinearConstraintDOF*				m_parentDof;	// parent degree of freedom
	std::vector<FELinearConstraintDOF*>	m_childDof;		// list of child dofs
	double			m_off;			// offset value

	DECLARE_PARAM_LIST();
};
