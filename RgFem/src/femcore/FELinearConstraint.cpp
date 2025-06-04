#include "FELinearConstraint.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "basicio/DumpStream.h"
#include "logger/log.h"

DEFINE_META_CLASS(FELinearConstraint, FEBoundaryCondition, "");
DEFINE_META_CLASS(FELinearConstraintDOF, FEObjectBase, "");

BEGIN_PARAM_DEFINE(FELinearConstraintDOF, FEObjectBase)
	ADD_PARAMETER(dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(node, "node");
	ADD_PARAMETER(val, "value");
END_PARAM_DEFINE();

FELinearConstraintDOF::FELinearConstraintDOF() : FEObjectBase()
{
	node = dof = -1;
	val = 1.0;
}

//=============================================================================
BEGIN_PARAM_DEFINE(FELinearConstraint, FEBoundaryCondition)
	ADD_PARAMETER(m_parentDof->dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(m_parentDof->node, "node");
	ADD_PARAMETER(m_off, "offset");

	//ADD_PROPERTY(m_childDof, "child_dof");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint() : FEBoundaryCondition()
{
	m_parentDof = new FELinearConstraintDOF();
	m_parentDof->GetParameterList(); // we need to call this to make sure that the parameter list is created.
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::~FELinearConstraint()
{
	Clear();
}

//-----------------------------------------------------------------------------
// return offset
double FELinearConstraint::GetOffset() const
{
	return m_off;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Clear()
{
	if (m_parentDof) delete m_parentDof; m_parentDof = nullptr;
	for (size_t i = 0; i < m_childDof.size(); ++i) delete m_childDof[i];
	m_childDof.clear();
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(const FELinearConstraint& LC) : FEBoundaryCondition()
{
	m_parentDof = nullptr;
	CopyFrom(&(const_cast<FELinearConstraint&>(LC)));
}

//-----------------------------------------------------------------------------
void FELinearConstraint::CopyFrom(FEBoundaryCondition* pbc)
{
	FELinearConstraint& LC = dynamic_cast<FELinearConstraint&>(*pbc);

	Clear();
	if (LC.m_parentDof)
	{
		m_parentDof = new FELinearConstraintDOF();
		m_parentDof->GetParameterList(); // NOTE: we need to call this to make sure that the parameter list is created.
		m_parentDof->node = LC.m_parentDof->node;
		m_parentDof->dof  = LC.m_parentDof->dof;
		m_parentDof->val  = LC.m_parentDof->val;
	}
	m_off = LC.m_off;
	int n = (int)LC.m_childDof.size();
	std::vector<FELinearConstraintDOF*>::const_iterator it = LC.m_childDof.begin();
	for (int i = 0; i < n; ++i, ++it)
	{
		FELinearConstraintDOF* d = new FELinearConstraintDOF();
		d->GetParameterList(); // NOTE: we need to call this to make sure that the parameter list is created.
		d->node = (*it)->node;
		d->dof  = (*it)->dof;
		d->val  = (*it)->val;
		m_childDof.push_back(d);
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentDof(int dof, int node)
{
	if (m_parentDof == nullptr) {
		m_parentDof = new FELinearConstraintDOF();
		m_parentDof->GetParameterList(); // NOTE: we need to call this to make sure that the parameter list is created.
	}
	m_parentDof->dof = dof;
	m_parentDof->node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentNode(int node)
{
	if (m_parentDof == nullptr) {
		m_parentDof = new FELinearConstraintDOF();
		m_parentDof->GetParameterList(); // NOTE: we need to call this to make sure that the parameter list is created.
	}
	m_parentDof->node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetParentDof(int dof)
{
	if (m_parentDof == nullptr) {
		m_parentDof = new FELinearConstraintDOF();
		m_parentDof->GetParameterList(); // NOTE: we need to call this to make sure that the parameter list is created.
	}
	m_parentDof->dof = dof;
}

//-----------------------------------------------------------------------------
// get the parent dof
int FELinearConstraint::GetParentDof() const
{
	return m_parentDof->dof;
}

//-----------------------------------------------------------------------------
int FELinearConstraint::GetParentNode() const
{
	return m_parentDof->node;
}

//-----------------------------------------------------------------------------
// get the child DOF
const FELinearConstraintDOF& FELinearConstraint::GetChildDof(int n) const
{
	return *m_childDof[n];
}

//-----------------------------------------------------------------------------
size_t FELinearConstraint::Size() const
{
	return m_childDof.size();
}

//-----------------------------------------------------------------------------
FELinearConstraint::dof_iterator FELinearConstraint::begin()
{
	return m_childDof.begin();
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddChildDof(int dof, int node, double v)
{
	FELinearConstraintDOF* d = new FELinearConstraintDOF();
	d->GetParameterList();	// we need to call this to make sure that the parameter list is created.
	d->dof = dof;
	d->node = node;
	d->val = v;
	m_childDof.push_back(d);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddChildDof(FELinearConstraintDOF* dof)
{
	dof->GetParameterList(); 	// we need to call this to make sure that the parameter list is created.
	m_childDof.push_back(dof);
}

//-----------------------------------------------------------------------------
// Initialization.
// Make sure the parent dof does not appear as a child dof
bool FELinearConstraint::Init()
{
	if (m_parentDof == nullptr) return false;
	if (m_parentDof->node < 0)
	{
		feLogError("Invalid node ID for parent node of linear constraint.");
		return false;
	}

	int n = (int)m_childDof.size();
	for (int i=0; i<n; ++i)
	{
		FELinearConstraintDOF& childNode = *m_childDof[i];
		if ((childNode.node == m_parentDof->node) && (childNode.dof == m_parentDof->dof)) return false;
		if (childNode.node < 0)
		{
			feLogError("Invalid node ID for child node  of linear constraint.");
			return false;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
// This is called during model activation (i.e. at the start of an analysis step)
// The parent dof is fixed in order to make sure that they are not assigned an equation number.
void FELinearConstraint::Activate()
{
	FEStepComponent::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	// we need the parent node to be fixed so that no equation is allocated
	FENode& node = mesh.Node(m_parentDof->node);
	node.setDofState(m_parentDof->dof, DOF_FIXED);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Deactivate()
{
	FEStepComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	FENode& node = mesh.Node(m_parentDof->node);
	node.setDofState(m_parentDof->dof, DOF_OPEN);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		m_parentDof->Serialize(ar);
		int n = (int)m_childDof.size();
		ar << n;
		std::vector<FELinearConstraintDOF*>::iterator it = m_childDof.begin();
		for (int i = 0; i < n; ++i, ++it) (*it)->Serialize(ar);
	}
	else
	{
		m_childDof.clear();
		if (m_parentDof == nullptr) {
			m_parentDof = new FELinearConstraintDOF(/*GetFEModel()*/);
			m_parentDof->GetParameterList(); // we need to call this to make sure that the parameter list is created.
		}
		m_parentDof->Serialize(ar);
		int n;
		ar >> n;
		for (int i=0; i<n; ++i)
		{
			FELinearConstraintDOF* dof = new FELinearConstraintDOF();
			dof->GetParameterList(); // we need to call this to make sure that the parameter list is created.
			dof->Serialize(ar);
			m_childDof.push_back(dof);
		}
	}
}
