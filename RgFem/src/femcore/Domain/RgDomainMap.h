#pragma once
#include "RgDomainList.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEDomainMap;

//-----------------------------------------------------------------------------
//! This class manages the domains in a model. It stores the domains in a list
//! and provides access via iterators. It also provides functionality for 
//! adding/removing domains. 
class FEM_EXPORT RgDomainMap : public RgDomainList
{
public:
	//! constructor
	RgDomainMap(FEModel* pfem);

	//! destructor
	~RgDomainMap();

	//! Create a new domain
	RgDomain* CreateDomain(const char* sztype);

	//! Add a domain to the list
	void AddDomain(RgDomain* pd);

public:
	// helper class for looping over all domains
	class Iterator
	{
	public:
		Iterator(RgDomainMap* pmap, int ndom = -1) : m_pmap(pmap), m_ndom(ndom) {}
		Iterator operator ++ () { m_ndom++; if (m_ndom >= m_pmap->Size()) m_ndom = -1; return (*this); }
		RgDomain* operator * () { return (m_ndom >= 0 ? (*m_pmap)[m_ndom] : 0); }
		bool operator != (const Iterator& it) { return (m_ndom != it.m_ndom); }
		Iterator begin() { return Iterator(m_pmap, 0); }
		Iterator end() { return Iterator(m_pmap, -1); }

	private:
		RgDomainMap*	m_pmap;
		int			m_ndom;
	};

	Iterator begin() { return Iterator(this, 0); }
	Iterator end() { return Iterator(this, -1); }

protected:
	void BuildDomainList();
};