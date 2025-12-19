#include "RgDomain.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "tools.h"
#include "log.h"

BEGIN_META_CLASS(RgDomain, 0)
	ADD_PROPERTY("name", &RgDomain::m_szname);
END_META_CLASS();

//-----------------------------------------------------------------------------
RgDomain::RgDomain(FE_Domain_Class dclass, FEModel* pfem) : m_dclass(dclass), m_pfem(pfem)
{
	m_pMesh = 0;
	m_pMat = 0;
	m_szname = "domain";
}

//-----------------------------------------------------------------------------
RgDomain::~RgDomain()
{
}

//-----------------------------------------------------------------------------
RgDomain::RgDomain(const RgDomain& d) : m_dclass(d.m_dclass), m_pfem(d.m_pfem)
{
	m_pMesh = d.m_pMesh;
	m_pMat = d.m_pMat;
	m_Node = d.m_Node;
	m_lnode = d.m_lnode;
	m_espec = d.m_espec;
	m_szname = d.m_szname;
}

//-----------------------------------------------------------------------------
RgDomain& RgDomain::operator = (const RgDomain& d)
{
	m_pMesh = d.m_pMesh;
	m_pMat = d.m_pMat;
	m_Node = d.m_Node;
	m_lnode = d.m_lnode;
	m_espec = d.m_espec;
	m_dclass = d.m_dclass;
	m_szname = d.m_szname;
	return (*this);
}

//-----------------------------------------------------------------------------
void RgDomain::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FECoreBase::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_lnode;
	}
	else
	{
		ar >> m_lnode;
	}
}

//-----------------------------------------------------------------------------
bool RgDomain::Init()
{
	// get the mesh
	m_pMesh = m_pfem->GetMesh();

	// allocate node list
	int NN = (int)m_lnode.size();
	m_Node.resize(NN);
	for (int i = 0; i < NN; ++i)
	{
		int nid = m_lnode[i];
		m_Node[i] = &m_pMesh->Node(nid);
	}

	return true;
}

//-----------------------------------------------------------------------------
void RgDomain::Reset()
{
}

//-----------------------------------------------------------------------------
void RgDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
}

//-----------------------------------------------------------------------------
void RgDomain::Activate()
{
}

//-----------------------------------------------------------------------------
void RgDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = pmat;
}

//-----------------------------------------------------------------------------
void RgDomain::ForEachElement(std::function<void(FEElement& el)> f)
{
}

//-----------------------------------------------------------------------------
int RgDomain::GetTotalDofs()
{
	return 0;
}