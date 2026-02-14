#include "RgElementSet.h"
#include "elements/RgElement/RgElement.h"
#include "femcore/FEMesh.h"
#include "femcore/Domain/RgDomain.h"
#include "basicio/DumpStream.h"
#include "femcore/FEModel.h"

//-----------------------------------------------------------------------------
RgElementSet::RgElementSet(FEModel* fem) : FEItemList(&fem->GetMesh())
{
	m_minID = -1;
	m_maxID = -1;
}

RgElementSet::RgElementSet(FEMesh* mesh) : FEItemList(mesh)
{
	m_minID = -1;
	m_maxID = -1;
}

//-----------------------------------------------------------------------------
void RgElementSet::Create(const std::vector<int>& elemList)
{
	//m_dom->Clear();
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void RgElementSet::Create(RgDomain* dom, const std::vector<int>& elemList)
{
	//m_dom.();
	m_dom->Add(dom);
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void RgElementSet::CopyFrom(RgElementSet& eset)
{
    /*SetName(eset.GetName());

    m_Elem = eset.m_Elem;

    FEMesh* mesh = GetMesh(); assert(mesh);

    m_dom.Clear();
    RgDomainList& dl = eset.GetDomainList();
    for (int i = 0; i < dl.Domains(); ++i)
    {
        RgDomain* di = dl.GetDomain(i);
        RgDomain* newdi = mesh->FindDomain(di->GetName()); assert(newdi);
        m_dom.AddDomain(newdi);
    }*/

	BuildLUT();
}

//-----------------------------------------------------------------------------
void RgElementSet::Create(RgDomain* dom)
{
	//m_dom.Clear();
	m_dom->Add(dom);
	SetMesh(dom->GetMesh());

	int NE = dom->Elements();
	m_Elem.resize(NE, -1);
	for (int i = 0; i < NE; ++i)
	{
		RgElement& el = dom->ElementRef(i);
		m_Elem[i] = el.getId();
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
// add another element set
void RgElementSet::Add(const RgElementSet& set)
{
	// add the domain list
	//m_dom.AddDomainList(set.GetDomainList());

	// add the elements
	m_Elem.insert(m_Elem.end(), set.m_Elem.begin(), set.m_Elem.end());

	BuildLUT();
}

//-----------------------------------------------------------------------------
void RgElementSet::Create(RgDomainList& domList)
{
	int NT = 0;
	//m_dom.Clear();
	for (int n = 0; n < domList.Size(); ++n)
	{
		RgDomain* dom = domList.Get(n);
		m_dom->Add(dom);
		NT += dom->Elements();
	}

	m_Elem.resize(NT, -1);
	NT = 0;
	for (int n = 0; n < domList.Size(); ++n)
	{
		RgDomain* dom = domList.Get(n);
		int NE = dom->Elements();
		for (int i = 0; i < NE; ++i)
		{
			RgElement& el = dom->ElementRef(i);
			m_Elem[NT + i] = el.getId();
		}
		NT += NE;
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void RgElementSet::BuildLUT()
{
	FEMesh* mesh = GetMesh();
	int N = (int)m_Elem.size();
	m_minID = m_maxID = -1;
	for (int i = 0; i < N; ++i)
	{
		RgElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->getId();

		if ((id < m_minID) || (m_minID == -1)) m_minID = id;
		if ((id > m_maxID) || (m_maxID == -1)) m_maxID = id;
	}

	int lutSize = m_maxID - m_minID + 1;
	m_LUT.resize(lutSize, -1);
	for (int i = 0; i < N; ++i)
	{
		RgElement* pe = mesh->FindElementFromID(m_Elem[i]);
		int id = pe->getId() - m_minID;
		m_LUT[id] = i;
	}
}

//-----------------------------------------------------------------------------
RgElement& RgElementSet::Element(int i)
{
	FEMesh* mesh = GetMesh();
	return *mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
const RgElement& RgElementSet::Element(int i) const
{
	FEMesh* mesh = GetMesh();
	return *mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
void RgElementSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Elem;
	ar & m_LUT;
	ar & m_minID & m_maxID;
}

void RgElementSet::SaveClass(DumpStream& ar, RgElementSet* p) {}
RgElementSet* RgElementSet::LoadClass(DumpStream& ar, RgElementSet* p)
{
	p = new RgElementSet(&ar.GetFEModel());
	return p;
}

//-----------------------------------------------------------------------------
// create node list from this element set
FENodeList RgElementSet::GetNodeList() const
{
	FEMesh* mesh = GetMesh();
	FENodeList set(mesh);
	std::vector<int> tag(mesh->Nodes(), 0);
	for (int i = 0; i<Elements(); ++i)
	{
		const RgElement& el = Element(i);
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