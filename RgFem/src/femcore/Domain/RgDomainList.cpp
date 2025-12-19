#include "RgDomainList.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
RgDomainList::RgDomainList(FEModel* pfem) : m_pfem(pfem)
{
}

//-----------------------------------------------------------------------------
RgDomainList::~RgDomainList()
{
	for (int i = 0; i < (int)m_Dom.size(); ++i) delete m_Dom[i];
	m_Dom.clear();
}

//-----------------------------------------------------------------------------
RgDomainList::RgDomainList(const RgDomainList& dl)
{
	m_pfem = dl.m_pfem;
	m_Dom = dl.m_Dom;
}

//-----------------------------------------------------------------------------
RgDomainList& RgDomainList::operator = (const RgDomainList& dl)
{
	m_pfem = dl.m_pfem;
	m_Dom = dl.m_Dom;
	return (*this);
}

//-----------------------------------------------------------------------------
void RgDomainList::Add(RgDomain* pd)
{
	m_Dom.push_back(pd);
}

//-----------------------------------------------------------------------------
void RgDomainList::Insert(int n, RgDomain* pd)
{
	m_Dom.insert(m_Dom.begin() + n, pd);
}

//-----------------------------------------------------------------------------
void RgDomainList::Remove(RgDomain* pd)
{
	for (int i = 0; i < (int)m_Dom.size(); ++i)
	{
		if (m_Dom[i] == pd)
		{
			m_Dom.erase(m_Dom.begin() + i);
			break;
		}
	}
}

//-----------------------------------------------------------------------------
bool RgDomainList::IsInList(RgDomain* pd)
{
	for (int i = 0; i < (int)m_Dom.size(); ++i)
	{
		if (m_Dom[i] == pd) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
void RgDomainList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		int ND = Size();
		ar << ND;
		for (int i = 0; i < ND; ++i) m_Dom[i]->Serialize(ar);
	}
	else
	{
		int ND;
		ar >> ND;
		for (int i = 0; i < ND; ++i) m_Dom[i]->Serialize(ar);
	}
}