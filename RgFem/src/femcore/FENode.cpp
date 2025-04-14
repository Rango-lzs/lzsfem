#include "FENode.h"
#include "basicio/DumpStream.h"

//=============================================================================
// FENode
//-----------------------------------------------------------------------------
FENode::FENode()
{
	// set the default state
	m_nstate = 0;

	// rigid body data
	m_rid = -1;

	// default ID
	mId = -1;
}

//-----------------------------------------------------------------------------
void FENode::SetDOFS(int n)
{
	// initialize dof stuff
	m_dofs.assign(n, -1);
	m_BC.assign(n, 0);
	m_val_t.assign(n, 0.0);
	m_val_p.assign(n, 0.0);
	m_Fr.assign(n, 0.0);
}

//-----------------------------------------------------------------------------
FENode::FENode(const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_d0 = n.m_d0;
    m_dt = n.m_dt;
    m_dp = n.m_dp;

	mId = n.mId;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_dofs = n.m_dofs;
	m_BC = n.m_BC;
	m_val_t = n.m_val_t;
	m_val_p = n.m_val_p;
	m_Fr = n.m_Fr;
}

//-----------------------------------------------------------------------------
FENode& FENode::operator = (const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_d0 = n.m_d0;
    m_dt = n.m_dt;
    m_dp = n.m_dp;

	mId = n.mId;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_dofs = n.m_dofs;
	m_BC = n.m_BC;
	m_val_t = n.m_val_t;
	m_val_p = n.m_val_p;
	m_Fr = n.m_Fr;

	return (*this);
}

//-----------------------------------------------------------------------------
// Serialize
void FENode::Serialize(DumpStream& ar)
{
	ar & mId;
	ar & m_rt & m_at;
	ar & m_rp & m_vp & m_ap;
	ar & m_Fr;
	ar & m_val_t & m_val_p;
    ar & m_dt & m_dp;
	if (ar.IsShallow() == false)
	{
		ar & m_nstate;
		ar & mId;
		ar & m_BC;
		ar & m_r0;
		ar & m_rid;
		ar & m_d0;
	}
}

//-----------------------------------------------------------------------------
//! Update nodal values, which copies the current values to the previous array
void FENode::UpdateValues()
{
	m_val_p = m_val_t;
}
