#include "FEDomainList.h"
#include "basicio/DumpStream.h"
#include "FEDomain.h"
#include "femcore/FEMesh.h"
#include <assert.h>
using namespace std;

FEDomainList::FEDomainList()
{

}

FEDomainList::FEDomainList(FEDomainList& domList)
{
	m_dom = domList.m_dom;
}

//! Clear the domain list
void FEDomainList::Clear()
{
	m_dom.clear();
}

void FEDomainList::AddDomain(FEDomain* dom)
{
	// see if this domain is already a member of this list
	if (IsMember(dom))
	{
//		assert(false);
		return;
	}

	// it's not, so let's add it
	m_dom.push_back(dom);
}

//! Add a domain list
void FEDomainList::AddDomainList(const FEDomainList& domList)
{
	for (int i = 0; i < domList.Domains(); ++i)
	{
		FEDomain* d = const_cast<FEDomain*>(domList.GetDomain(i));
		AddDomain(d);
	}
}

bool FEDomainList::IsMember(const FEDomain* dom) const
{
	// loop over all the domains
	for (size_t i = 0; i < m_dom.size(); ++i)
	{
		if (m_dom[i] == dom)
		{
			// found it!
			return true;
		}
	}

	// better luck next time!
	return false;
}

//! serialization
void FEDomainList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_dom;
}
