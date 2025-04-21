#pragma once

#include "femcore/FEItemList.h"
#include "elements/RgElement.h"
#include "femcore/Domain/FEDomainList.h"
#include "femcore/FENodeList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;
class DumpStream;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FEM_EXPORT FEElementSet : public FEItemList
{
public:
	//! constructor
	FEElementSet(FEModel* fem);	// TODO: remove!
	FEElementSet(FEMesh* fem);

	// Create the element set
	void Create(const std::vector<int>& elemList);
	void Create(FEDomain* dom, const std::vector<int>& elemList);

	void CopyFrom(FEElementSet& eset);

	// add another element set
	void Add(const FEElementSet& set);

	// Create the element set from a domain
	void Create(FEDomain* dom);

	// Create the element set from a domain
	void Create(FEDomainList& dom);

	// Return number of elements in the set
	int Elements() const { return (int)m_Elem.size(); }

	int operator [] (int i) const { return m_Elem[i]; }

	// return the local index of an element into the element set
	// returns -1 if the element is not part of element set
	int GetLocalIndex(const FEElement& el) const;
	bool Contains(const FEElement& el) const;

	// Get the element ID list
	const std::vector<int>& GetElementIDList() const { return m_Elem; }

	// get the domain list that generated the element set
	FEDomainList& GetDomainList() { return m_dom; }
	const FEDomainList& GetDomainList() const { return m_dom; }

	// Get an element
	FEElement& Element(int i);
	const FEElement& Element(int i) const;

	// create node list from this element set
	FENodeList GetNodeList() const;

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FEElementSet* p);
	static FEElementSet* LoadClass(DumpStream& ar, FEElementSet* p);

private:
	// Build the lookup table
	void BuildLUT();

protected:
	std::vector<int>	m_Elem;		//!< list of elements' global ID

	FEDomainList		m_dom;	//!< domain list that generated the element set

	// used for fast lookup in GetLocalIndex
	std::vector<int>	m_LUT;
	int					m_minID, m_maxID;
};

inline int FEElementSet::GetLocalIndex(const FEElement& el) const
{
    int eid = el.getID();
	if ((eid < m_minID) || (eid > m_maxID)) return -1;
    else
        return m_LUT[el.getID() - m_minID];
}

inline bool FEElementSet::Contains(const FEElement& el) const
{
	return (GetLocalIndex(el) != -1);
}
