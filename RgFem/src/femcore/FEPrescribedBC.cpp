#include "FEPrescribedBC.h"
#include "FESurface.h"
#include "FENode.h"

DEFINE_META_CLASS(FEPrescribedNodeSet, FENodalBC, "");

//=============================================================================
BEGIN_PARAM_DEFINE(FEPrescribedNodeSet, FENodalBC)
	ADD_PARAMETER(m_brelative, "relative");
END_PARAM_DEFINE();

FEPrescribedNodeSet::FEPrescribedNodeSet() : FENodalBC()
{
	m_brelative = false;
}

//-----------------------------------------------------------------------------
// set the relative flag
void FEPrescribedNodeSet::SetRelativeFlag(bool br)
{
	m_brelative = br;
}

void FEPrescribedNodeSet::Activate()
{
	FENodalBC::Activate();

	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	if (m_brelative) m_rval.assign(N * dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// set the dofs to prescribed
		for (size_t j = 0; j < dofs; ++j)
		{
			node.setDofState(m_dof[j], DOF_PRESCRIBED);
			if (m_brelative)
			{
				m_rval[i * dofs + j] = node.get(m_dof[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedNodeSet::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// set the dof to open
		for (int j = 0; j < dofs; ++j)
		{
			node.setDofState(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedNodeSet::PrepStep(std::vector<double>& ui, bool brel)
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			int I = -node.getDofs()[m_dof[j]] - 2;
			if (I >= 0) ui[I] = (brel ? uj - node.get(m_dof[j]) : uj);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedNodeSet::Serialize(DumpStream& ar)
{
	FENodalBC::Serialize(ar);
	ar & m_rval;
}

//-----------------------------------------------------------------------------
// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedNodeSet::Update()
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			node.set(m_dof[j], uj);
		}
	}
}

//-----------------------------------------------------------------------------
// This is called during contact update and should be used to enforce the
// nodal degrees of freedoms
void FEPrescribedNodeSet::Repair()
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			if (node.getDofs()[m_dof[j]] >= 0) {
                node.setDofIdx(m_dof[j], -node.getDofs()[m_dof[j]] - 2);
				double uj = val[j];
				if (m_brelative)
				{
					uj += m_rval[i * dofs + j];
				}

				node.set(m_dof[j], uj);
			}
		}
	}
}

//=============================================================================

BEGIN_PARAM_DEFINE(FEPrescribedSurface, FESurfaceBC)
END_PARAM_DEFINE();

FEPrescribedSurface::FEPrescribedSurface() : FESurfaceBC()
{
	m_brelative = false;
}

void FEPrescribedSurface::Activate()
{
	FESurfaceBC::Activate();

	FESurface* surface = GetSurface();
	if (surface == nullptr) return;

	m_nodeList = surface->GetNodeList();
	
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	if (m_brelative) m_rval.assign(N * dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// set the dofs to prescribed
		for (size_t j = 0; j < dofs; ++j)
		{
			node.setDofState(m_dof[j], DOF_PRESCRIBED);

			if (m_brelative)
			{
				m_rval[i * dofs + j] = node.get(m_dof[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedSurface::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// set the dof to open
		for (int j = 0; j < dofs; ++j)
		{
			node.setDofState(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedSurface::PrepStep(std::vector<double>& ui, bool brel)
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			int I = -node.getDofs()[m_dof[j]] - 2;
			if (I >= 0) ui[I] = (brel ? uj - node.get(m_dof[j]) : uj);
		}
	}
}

//-----------------------------------------------------------------------------
// set the relative flag
void FEPrescribedSurface::SetRelativeFlag(bool br)
{
	m_brelative = br;
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedSurface::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	ar & m_rval;
	if (ar.IsShallow() == false) ar & m_nodeList;
}

//-----------------------------------------------------------------------------
// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedSurface::Update()
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			node.set(m_dof[j], uj);
		}
	}
}

//-----------------------------------------------------------------------------
// This is called during contact update and should be used to enforce the
// nodal degrees of freedoms
void FEPrescribedSurface::Repair()
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			if (node.getDofs()[m_dof[j]] >= 0) {
				node.setDofIdx(m_dof[j],  -node.getDofs()[m_dof[j]] - 2);
				double uj = val[j];
				if (m_brelative)
				{
					uj += m_rval[i * dofs + j];
				}

				node.set(m_dof[j], uj);
			}
		}
	}
}
