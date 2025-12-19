#pragma once
#include "RgDomain.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;

//-----------------------------------------------------------------------------
//! A class that manages a list of domains
class FEM_EXPORT RgDomainList
{
public:
	//! constructor
	RgDomainList(FEModel* pfem);

	//! destructor
	~RgDomainList();

	//! Copy constructor
	RgDomainList(const RgDomainList& dl);

	//! Assignment operator
	RgDomainList& operator = (const RgDomainList& dl);

	//! add a domain to the list
	void Add(RgDomain* pd);

	//! insert a domain at a specific location
	void Insert(int n, RgDomain* pd);

	//! remove a domain
	void Remove(RgDomain* pd);

	//! return nr of domains
	int Size() const { return (int)m_Dom.size(); }

	//! return a domain
	RgDomain* Get(int i) { return m_Dom[i]; }

	//! return a domain
	const RgDomain* Get(int i) const { return m_Dom[i]; }

	//! operator overloading
	RgDomain* operator [] (int i) { return m_Dom[i]; }

	//! const operator overloading
	const RgDomain* operator [] (int i) const { return m_Dom[i]; }

	//! see if a domain is in the list
	bool IsInList(RgDomain* pd);

public:
	//! serialization
	void Serialize(DumpStream& ar);

protected:
	std::vector<RgDomain*>	m_Dom;
	FEModel*		m_pfem;
};