#include "FEElementSet.h"
#include "elements/RgElement.h"
#include "femcore/FEMesh.h"
#include "femcore/FEDomain.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEModel* fem) : FEItemList(fem)
{
	m_minID = -1;
	m_maxID = -1;
}

FEElementSet::FEElementSet(FEMesh* mesh) : FEItemList(mesh)
{
	m_minID = -1;
	m_maxID = -1;
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(const std::vector<int>& elemList)
{
	m_dom.Clear();
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomain* dom, const std::vector<int>& elemList)
{
	m_dom.Clear();
	m_dom.AddDomain(dom);
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::CopyFrom(FEElementSet& eset)
{
	SetName(eset.GetName());

	m_Elem = eset.m_Elem;

	FEMesh* mesh = GetMesh(); assert(mesh);

	m_dom.Clear();
	FEDomainList& dl = eset.GetDomainList();
	for (int i = 0; i < dl.Domains(); ++i)
	{
		FEDomain* di = dl.GetDomain(i);
		FEDomain* newdi = mesh->FindDomain(di->GetName()); assert(newdi);
		m_dom.AddDomain(newdi);
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomain* dom)
{
	m_dom.Clear();
	m_dom.AddDomain(dom);
	SetMesh(dom->GetMesh());

	int NE = dom->Elements();
	m_Elem.resize(NE, -1);
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom->ElementRef(i);
		m_Elem[i] = el.getID();
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
// add another element set
void FEElementSet::Add(const FEElementSet& set)
{
	// add the domain list
	m_dom.AddDomainList(set.GetDomainList());

	// add the elements
	m_Elem.insert(m_Elem.end(), set.m_Elem.begin(), set.m_Elem.end());

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomainList& domList)
{
	int NT = 0;
	m_dom.Clear();
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		m_dom.AddDomain(dom);
		NT += dom->Elements();
	}

	m_Elem.resize(NT, -1);
	NT = 0;
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		int NE = dom->Elements();
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = dom->ElementRef(i);
			m_Elem[NT + i] = el.getID();
		}
		NT += NE;
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::BuildLUT()
{
	FEMesh* mesh = GetMesh();
	int N = (int)m_Elem.size();
	m_minID = m_maxID = -1;
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->getID();

		if ((id < m_minID) || (m_minID == -1)) m_minID = id;
		if ((id > m_maxID) || (m_maxID == -1)) m_maxID = id;
	}

	int lutSize = m_maxID - m_minID + 1;
	m_LUT.resize(lutSize, -1);
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->getID() - m_minID;
		m_LUT[id] = i;
	}
}

//-----------------------------------------------------------------------------
FEElement& FEElementSet::Element(int i)
{
	FEMesh* mesh = GetMesh();
	return *mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
const FEElement& FEElementSet::Element(int i) const
{
	FEMesh* mesh = GetMesh();
	return *mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
void FEElementSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Elem;
	ar & m_LUT;
	ar & m_minID & m_maxID;
}

void FEElementSet::SaveClass(DumpStream& ar, FEElementSet* p) {}
FEElementSet* FEElementSet::LoadClass(DumpStream& ar, FEElementSet* p)
{
	p = new FEElementSet(&ar.GetFEModel());
	return p;
}

//-----------------------------------------------------------------------------
// create node list from this element set
FENodeList FEElementSet::GetNodeList() const
{
	FEMesh* mesh = GetMesh();
	FENodeList set(mesh);
	std::vector<int> tag(mesh->Nodes(), 0);
	for (int i = 0; i<Elements(); ++i)
	{
		const FEElement& el = Element(i);
		int ne = el.NodeSize();
		for (int j = 0; j<ne; ++j)
		{
			/*if (tag[el.m_node[j]] == 0)
			{
				set.Add(el.m_node[j]);
				tag[el.m_node[j]] = 1;
			}*/
		}
	}
	return set;
}
