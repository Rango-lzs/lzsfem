#include "FEMeshPartition.h"
#include "materials/FEMaterial.h"
//#include "FEDataExport.h"
#include "FEMesh.h"
#include "DOFS.h"
#include <string.h>
#include "FEModel.h"
#include "basicio/DumpStream.h"

 DEFINE_META_CLASS(FEMeshPartition, FEObjectBase, "");

//-----------------------------------------------------------------------------
FEMeshPartition::FEMeshPartition(int nclass, FEModel* fem) : FEObjectBase(), m_nclass(nclass)
{
	m_pMesh = nullptr;
	if (fem) m_pMesh = &fem->GetMesh();
	m_bactive = true;
}

//-----------------------------------------------------------------------------
FEMeshPartition::~FEMeshPartition()
{
	// delete all data export classes
	if (m_Data.empty() == false)
	{
		size_t ND = m_Data.size();
		for (size_t i = 0; i<ND; ++i) delete m_Data[i];
		m_Data.clear();
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::AddDataExport(FEDataExport* pd)
{
	if (pd) m_Data.push_back(pd);
}

//-----------------------------------------------------------------------------
FEElement* FEMeshPartition::FindElementFromID(int nid)
{
	for (int i = 0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.getId() == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------
void FEMeshPartition::Serialize(DumpStream& ar)
{
	FEObjectBase::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Node;
	ar & m_nclass;
	ar & m_bactive;
//	ar & m_Data;
}

//-----------------------------------------------------------------------------
//! return a specific node
FENode& FEMeshPartition::Node(int i)
{
	return m_pMesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
//! return a specific node
const FENode& FEMeshPartition::Node(int i) const
{
	return m_pMesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
void FEMeshPartition::CopyFrom(FEMeshPartition* pd)
{
	m_Node = pd->m_Node;
	SetName(pd->GetName());
}

//-----------------------------------------------------------------------------
bool FEMeshPartition::Init()
{
	// base class first
	if (FEObjectBase::Init() == false) return false;

	// make sure that there are elements in this domain
	if (Elements() == 0) return false;

	// get the mesh to which this domain belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	int NN = mesh.Nodes();
	std::vector<int> tag; tag.assign(NN, -1);

	// let's find all nodes the domain needs
	int nn = 0;
	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& el = ElementRef(i);
		int ne = el.NodeSize();
		for (int j = 0; j<ne; ++j)
		{
			// get the global node number
			int m = el.m_node[j];

			// create a local node number
			if (tag[m] == -1) tag[m] = nn++;

			// set the local node number
			el.m_loc_node[j] = tag[m];
		}
	}

	// allocate node index table
	m_Node.assign(nn, -1);

	// fill the node index table
	for (int i = 0; i<NN; ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}

#ifdef _DEBUG
	// make sure all nodes are assigned a local index
	for (int i = 0; i<nn; ++i)
	{
		assert(m_Node[i] >= 0);
	}
#endif

	return true;
}


//-----------------------------------------------------------------------------
void FEMeshPartition::ForEachMaterialPoint(std::function<void(FEMaterialPoint& mp)> f)
{
	int NE = Elements();
#pragma omp parallel for shared(f)
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = ElementRef(i);
		int nint = el.GaussPointSize();
		for (int n = 0; n < nint; ++n) f(*el.GetMaterialPoint(n));
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::ForEachElement(std::function<void(FEElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}
