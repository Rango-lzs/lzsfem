#pragma once

#include "femcore/FEItemList.h"
#include "elements/RgElement/RgElement.h"
#include "femcore/Domain/RgDomainList.h"
#include "femcore/FENodeList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;
class DumpStream;
class FEModel;
class RgDomainList;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FEM_EXPORT RgElementSet : public FEItemList
{
public:
	//! constructor
	RgElementSet(FEModel* fem);	// TODO: remove!
	RgElementSet(FEMesh* fem);

	// Create the element set
	void Create(const std::vector<int>& elemList);
	void Create(RgDomain* dom, const std::vector<int>& elemList);

	void CopyFrom(RgElementSet& eset);

	// add another element set
	void Add(const RgElementSet& set);

	// Create the element set from a domain
	void Create(RgDomain* dom);

	// Create the element set from a domain
	void Create(RgDomainList& dom);

	// Return number of elements in the set
	int Elements() const { return (int)m_Elem.size(); }

	int operator [] (int i) const { return m_Elem[i]; }

	// return the local index of an element into the element set
	// returns -1 if the element is not part of element set
	int GetLocalIndex(const RgElement& el) const;
	bool Contains(const RgElement& el) const;

	// Get the element ID list
	const std::vector<int>& GetElementIDList() const { return m_Elem; }

	// get the domain list that generated the element set
	RgDomainList* GetDomainList() { return m_dom; }
	const RgDomainList* GetDomainList() const { return m_dom; }

	// Get an element
	RgElement& Element(int i);
	const RgElement& Element(int i) const;

	// create node list from this element set
	FENodeList GetNodeList() const;

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, RgElementSet* p);
	static RgElementSet* LoadClass(DumpStream& ar, RgElementSet* p);

private:
	// Build the lookup table
	void BuildLUT();

protected:
	std::vector<int>	m_Elem;		//!< list of elements' global ID

	RgDomainList*		m_dom;	//!< domain list that generated the element set

	// used for fast lookup in GetLocalIndex
	std::vector<int>	m_LUT;
	int					m_minID, m_maxID;
};

inline int RgElementSet::GetLocalIndex(const RgElement& el) const
{
    int eid = el.getId();
	if ((eid < m_minID) || (eid > m_maxID)) return -1;
    else
        return m_LUT[el.getId() - m_minID];
}

inline bool RgElementSet::Contains(const RgElement& el) const
{
	return (GetLocalIndex(el) != -1);
}