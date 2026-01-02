#include "RgDomain.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENode.h"
#include "femcore/tools.h"
#include "logger/log.h"
#include "../Matrix/FEGlobalMatrix.h"


//-----------------------------------------------------------------------------
RgDomain::RgDomain(FEModel* pfem) : m_pfem(pfem)
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
RgDomain::RgDomain(const RgDomain& d) : m_pfem(d.m_pfem)
{
	m_pMesh = d.m_pMesh;
	m_pMat = d.m_pMat;
	m_Node = d.m_Node;
	m_lnode = d.m_lnode;
	m_szname = d.m_szname;
}

//-----------------------------------------------------------------------------
RgDomain& RgDomain::operator = (const RgDomain& d)
{
	m_pMesh = d.m_pMesh;
	m_pMat = d.m_pMat;
	m_Node = d.m_Node;
	m_lnode = d.m_lnode;
	m_szname = d.m_szname;
	return (*this);
}

//-----------------------------------------------------------------------------
void RgDomain::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEObjectBase::Serialize(ar);

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
	m_pMesh = &m_pfem->GetMesh();

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


void RgDomain::BuildMatrixProfile(FEGlobalMatrix& M)
{
    std::vector<int> elm;
    const int NE = Elements();
    for (int j = 0; j < NE; ++j)
    {
        RgElement& el = ElementRef(j);
        UnpackLM(el, elm);
        M.build_add(elm);
    }
}

//-----------------------------------------------------------------------------
void RgDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = pmat;
}

//-----------------------------------------------------------------------------
void RgDomain::ForEachElement(std::function<void(RgElement& el)> f)
{
    const int n = Elements();
    for (int i = 0; i < n; ++i)
    {
        f(ElementRef(i));
    }
}

//-----------------------------------------------------------------------------
int RgDomain::GetTotalDofs()
{
    int dofs = 0;
    const int n = Nodes();
    for (int i = 0; i < n; ++i)
    {
        FENode& node = Node(i);
        dofs += node.dofSize();
    }
    return dofs;
}