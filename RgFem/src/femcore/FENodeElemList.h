#pragma once
#include "femcore/fem_export.h"
#include <vector>

class FESurface;
class FEMesh;
class RgElement;
class RgDomain;
class DumpStream;

//-----------------------------------------------------------------------------
//! The FENodeElemList class is a utility class that determines for each node 
//! to which element it belongs.

//! This class analyzes a mesh and finds for each node all elements that have
//! this node

class FEM_EXPORT FENodeElemList
{
public:
	FENodeElemList(){}
	virtual ~FENodeElemList(){}

	//! build the node-element list for a surface
	void Create(const FESurface& s);

	//! build the node-selement list for a mesh
	void Create(FEMesh& mesh);

	//! build the node-element list for a domain
	void Create(RgDomain& dom);

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar);

	//! Clear the list
	void Clear();

	int MaxValence();
	int Valence(int n) { return m_nval[n]; }
	RgElement** ElementList(int n) { return &m_eref[0] + m_pn[n]; }
	int* ElementIndexList(int n) { return &m_iref[0] + m_pn[n]; }

	int Size() { return (int) m_nval.size(); }

protected:
	std::vector<int>			m_nval;	// nodal valences
	std::vector<RgElement*>		m_eref;	// element pointers
	std::vector<int>			m_iref;	// element indices
	std::vector<int>			m_pn;	// start index into the eref array
};

//-----------------------------------------------------------------------------
//! Like the FEElemElemList, but can create multiple levels
class FENodeElemTree
{
public:
	FENodeElemTree() {}
	virtual ~FENodeElemTree() {}

	void Create(FESurface* ps, int k = 0);

	int Valence(int n) { return (int) m_nel[n].size(); }

	RgElement** ElementList(int n) { return &(m_nel[n][0]);}

	bool empty() { return m_nel.empty(); }

protected:
	std::vector< std::vector<RgElement*> >	m_nel;
};
