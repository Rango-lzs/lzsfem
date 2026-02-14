#include "FEMesh.h"

#include "FEException.h"
// #include "femcore/Domain/FEDiscreteDomain.h"
#include "femcore/Domain/RgTrussDomain.h"
// #include "femcore/Domain/FEShellDomain.h"
#include "elements/RgElemElemList.h"
#include "elements/RgElementList.h"
#include "FEDataArray.h"
#include "femcore/DOFS.h"
#include "femcore/Domain/RgSolidDomain.h"
#include "femcore/FESurface.h"
// #include "femcore/Domain/RgDomainMap.h"
#include "basicio/DumpStream.h"
#include "Domain/RgDomain.h"
#include "FENodeDataMap.h"
#include "FESurfaceMap.h"

#include <algorithm>

//-----------------------------------------------------------------------------
FEDataMap* CreateDataMap(int mapType)
{
    FEDataMap* map = nullptr;
    switch (mapType)
    {
        case FE_NODE_DATA_MAP:
            map = new FENodeDataMap;
            break;
        // case FE_DOMAIN_MAP   : map = new FEDomainMap; break;
        case FE_SURFACE_MAP:
            map = new FESurfaceMap;
            break;
        default:
            assert(false);
    }

    return map;
}

//=============================================================================
// FEMesh
//-----------------------------------------------------------------------------
FEMesh::FEMesh(FEModel* fem)
    : m_fem(fem)
{
    m_LUT = 0;
}

//-----------------------------------------------------------------------------
FEMesh::~FEMesh()
{
    Clear();
}

//-----------------------------------------------------------------------------
//! return number of nodes
int FEMesh::Nodes() const
{
    return (int)m_Node.size();
}

//-----------------------------------------------------------------------------
const std::vector<FENode>& FEMesh::AllNode() const
{
    return m_Node;
}

//-----------------------------------------------------------------------------
FEBoundingBox& FEMesh::GetBoundingBox()
{
    return m_box;
}

//-----------------------------------------------------------------------------
int FEMesh::NodeSets() const
{
    return (int)m_NodeSet.size();
}

//-----------------------------------------------------------------------------
int FEMesh::ElementSets() const
{
    return (int)m_ElemSet.size();
}

//-----------------------------------------------------------------------------
const std::vector<RgDomain*>& FEMesh::AllDomain() const
{
    return m_Domain;
}

//-----------------------------------------------------------------------------
int FEMesh::Surfaces() const
{
    return (int)m_Surf.size();
}

//-----------------------------------------------------------------------------
int FEMesh::Edges() const
{
    return (int)m_Edge.size();
}

//-----------------------------------------------------------------------------
int FEMesh::FacetSets() const
{
    return (int)m_FaceSet.size();
}

//-----------------------------------------------------------------------------
int FEMesh::SurfacePairs() const
{
    return (int)m_SurfPair.size();
}

//-----------------------------------------------------------------------------
void FEMesh::Serialize(DumpStream& ar)
{
    // clear the mesh if we are loading from an archive
    if ((ar.IsShallow() == false) && (ar.IsLoading()))
        Clear();

    // we don't want to store pointers to all the nodes
    // mostly for efficiency, so we tell the archive not to store the pointers
    ar.LockPointerTable();
    {
        // store the node list
        ar& m_Node;
    }
    ar.UnlockPointerTable();

    // stream domain data
    ar& m_Domain;

    // if this is a shallow archive, we're done
    if (ar.IsShallow())
        return;

    // serialize node sets
    ar& m_NodeSet;

    if (ar.IsSaving())
    {
        // write segment sets
        // Rango TODO:
        /*int ssets = SegmentSets();
        ar << ssets;
        for (int i=0; i<ssets; ++i)
        {
            FESegmentSet& sset = SegmentSet(i);
            sset.Serialize(ar);
        }*/

        // write element sets
        ar << m_ElemSet;

        // write facet sets
        ar << m_FaceSet;

        // write surface pairs
        int surfPairs = m_SurfPair.size();
        ar << surfPairs;
        for (int i = 0; i < m_SurfPair.size(); ++i)
        {
            FESurfacePair& sp = *m_SurfPair[i];
            ar << sp.GetName();
            ar << sp.GetPrimarySurface()->GetName();
            ar << sp.GetSecondarySurface()->GetName();
        }

        // write discrete sets
        /*int dsets = DiscreteSets();
        ar << dsets;
        for (int i=0; i<dsets; ++i)
        {
            FEDiscreteSet& dset = DiscreteSet(i);
            dset.Serialize(ar);
        }*/

        //// write data maps
        // int maps = DataMaps();
        // ar << maps;
        // for (int i = 0; i < maps; ++i)
        //{
        //	FEDataMap* map = GetDataMap(i);
        //	ar << map;
        // }
    }
    else
    {
        FEModel* fem = &ar.GetFEModel();

        // read segment sets
        int ssets = 0;
        ar >> ssets;
        /*for (int i=0; i<ssets; ++i)
        {
            FESegmentSet* sset = new FESegmentSet(fem);
            AddSegmentSet(sset);
            sset->Serialize(ar);
        }*/

        // read element sets
        ar >> m_ElemSet;

        // read facet sets
        ar >> m_FaceSet;

        // read surface pairs
        int surfPairs = 0;
        ar >> surfPairs;
        for (int i = 0; i < surfPairs; ++i)
        {
            FESurfacePair* sp = new FESurfacePair(this);
            std::string name;
            ar >> name;
            sp->SetName(name);

            ar >> name;
            FEFacetSet* ps = FindFacetSet(name);
            sp->SetPrimarySurface(ps);

            ar >> name;
            ps = FindFacetSet(name);
            sp->SetSecondarySurface(ps);

            AddSurfacePair(sp);
        }

        // read discrete sets
        int dsets = 0;
        ar >> dsets;
        for (int i = 0; i < dsets; ++i)
        {
            /*FEDiscreteSet* dset = new FEDiscreteSet(this);
            AddDiscreteSet(dset);
            dset->Serialize(ar);*/
        }

        //// write data maps
        // ClearDataMaps();
        // int maps = 0;
        // std::string mapName;
        // ar >> maps;
        // for (int i = 0; i < maps; ++i)
        //{
        //	FEDataMap* map = nullptr;
        //	ar >> map;
        //	AddDataMap(map);
        // }

        UpdateBox();
    }
}

//-----------------------------------------------------------------------------
void FEMesh::SaveClass(DumpStream& ar, FEMesh* p)
{
    // we should never get here
    assert(false);
}

//-----------------------------------------------------------------------------
FEMesh* FEMesh::LoadClass(DumpStream& ar, FEMesh* p)
{
    // we should never get here
    assert(false);
    return nullptr;
}

//-----------------------------------------------------------------------------
//  Allocates storage for mesh data.
//
void FEMesh::CreateNodes(int nodes)
{
    assert(nodes);
    m_Node.resize(nodes);

    // set the default node IDs
    for (int i = 0; i < nodes; ++i)
        Node(i).SetID(i + 1);

    m_NEL.Clear();
}

//-----------------------------------------------------------------------------
// Make more room for nodes
void FEMesh::AddNodes(int nodes)
{
    assert(nodes);
    int N0 = (int)m_Node.size();

    // get the ID of the last node
    // (It is assumed that nodes are sorted according their ID
    //  so the last node should have the highest ID)
    int n0 = 1;
    if (N0 > 0)
        n0 = m_Node[N0 - 1].GetID() + 1;

    m_Node.resize(N0 + nodes);
    for (int i = 0; i < nodes; ++i)
        m_Node[i + N0].SetID(n0 + i);
}

//-----------------------------------------------------------------------------
void FEMesh::SetDOFS(int n)
{
    int NN = Nodes();
    for (int i = 0; i < NN; ++i)
        m_Node[i].SetDOFS(n);
}

//-----------------------------------------------------------------------------
//! Return the total number elements
int FEMesh::Elements() const
{
    int N = 0;
    for (int i = 0; i < (int)m_Domain.size(); ++i)
    {
        N += m_Domain[i]->Elements();
    }
    return N;
}

//-----------------------------------------------------------------------------
//! Return the total number of elements of a specific domain type
int FEMesh::Elements(int ndom_type) const
{
    int N = 0;
    for (int i = 0; i < (int)m_Domain.size(); ++i)
    {
        RgDomain& dom = *m_Domain[i];
        if (dom.domType() == ndom_type)
            N += m_Domain[i]->Elements();
    }
    return N;
}

//-----------------------------------------------------------------------------
//! return reference to a node
FENode& FEMesh::Node(int i)
{
    return m_Node[i];
}
const FENode& FEMesh::Node(int i) const
{
    return m_Node[i];
}

//-----------------------------------------------------------------------------
//  Updates the bounding box of the mesh (using current coordinates)
//
void FEMesh::UpdateBox()
{
    if (Nodes() > 0)
    {
        m_box = FEBoundingBox(Node(0).m_rt);
        for (int i = 1; i < Nodes(); ++i)
        {
            m_box.add(Node(i).m_rt);
        }
    }
    else
        m_box = FEBoundingBox(Vector3d(0, 0, 0));
}

//-----------------------------------------------------------------------------
//  Counts the number of shell elements in the mesh
//
int FEMesh::RemoveIsolatedVertices()
{
    int i, j, k, N = Nodes(), n;

    // create a valence array
    std::vector<int> val;
    val.assign(N, 0);

    // count the nodal valences
    for (i = 0; i < (int)m_Domain.size(); ++i)
    {
        RgDomain& d = Domain(i);
        for (j = 0; j < d.Elements(); ++j)
        {
            RgElement& el = d.ElementRef(j);
            n = el.NodeSize();
            for (k = 0; k < n; ++k)
                ++val[el.getNodeId(k)];
        }
    }

    // See if there are any isolated nodes
    // Exclude them from the analysis
    int ni = 0;
    for (i = 0; i < N; ++i)
        if (val[i] == 0)
        {
            ++ni;
            FENode& node = Node(i);
            node.SetFlags(FENode::EXCLUDE);
        }

    return ni;
}

//-----------------------------------------------------------------------------
//! Does one-time initialization of the Mesh material point data.
void FEMesh::InitMaterialPoints()
{
    for (int i = 0; i < Domains(); ++i)
    {
        RgDomain& dom = Domain(i);
        // dom.InitMaterialPoints();
    }
}

//-----------------------------------------------------------------------------
void FEMesh::Clear()
{
    m_Node.clear();
    for (size_t i = 0; i < m_Domain.size(); ++i)
        delete m_Domain[i];

    // TODO: Surfaces are currently managed by the classes that use them so don't delete them
    //	for (size_t i=0; i<m_Surf.size   (); ++i) delete m_Surf   [i];

    for (size_t i = 0; i < m_NodeSet.size(); ++i)
        delete m_NodeSet[i];
    // for (size_t i=0; i<m_LineSet.size (); ++i) delete m_LineSet [i];
    for (size_t i = 0; i < m_ElemSet.size(); ++i)
        delete m_ElemSet[i];
    for (size_t i = 0; i < m_DiscSet.size(); ++i)
        delete m_DiscSet[i];
    for (size_t i = 0; i < m_FaceSet.size(); ++i)
        delete m_FaceSet[i];
    for (size_t i = 0; i < m_SurfPair.size(); ++i)
        delete m_SurfPair[i];

    m_Domain.clear();
    m_Surf.clear();
    //m_NodeSet.clear();
    // m_LineSet.clear();
    m_ElemSet.clear();
    m_DiscSet.clear();
    m_FaceSet.clear();
    m_SurfPair.clear();

    m_NEL.Clear();
    if (m_LUT)
        delete m_LUT;
    m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Reset the mesh data. Return nodes to their intial position, reset their
//! attributes and zero all element stresses.

void FEMesh::Reset()
{
    // reset nodal data
    for (int i = 0; i < Nodes(); ++i)
    {
        FENode& node = Node(i);

        node.m_rp = node.m_rt = node.m_r0;
        node.m_vp = Vector3d(0, 0, 0);
        node.m_ap = node.m_at = Vector3d(0, 0, 0);
        node.m_dp = node.m_dt = node.m_d0;

        // reset ID arrays
        int ndof = (int)node.dofSize();
        for (int i = 0; i < ndof; ++i)
        {
            node.set_inactive(i);  // Rango
            node.setDofState(i, DOF_OPEN);
            node.set(i, 0.0);
            node.set_load(i, 0.0);
        }
    }

    // update the mesh
    UpdateBox();

    // reset domain data
    for (int n = 0; n < (int)m_Domain.size(); ++n)
        m_Domain[n]->Reset();
}

//-----------------------------------------------------------------------------
//! This function calculates the (initial) volume of an element. In some case, the volume
//! may only be approximate.
double FEMesh::ElementVolume(RgElement& el)
{
    double V = 0;
    switch (el.elementType())
    {
        case FE_ELEM_SOLID:
        {
            RgSolidDomain* dom = dynamic_cast<RgSolidDomain*>(el.getDomain());
            assert(dom);
            // if (dom) V = dom->Volume(static_cast<RgSolidElement&>(el));
        }
        break;
        case FE_ELEM_SHELL:
        {
            /*FEShellDomain* dom = dynamic_cast<FEShellDomain*>(el.GetMeshPartition()); assert(dom);
            if (dom) V = dom->Volume(static_cast<FEShellElement&>(el));*/
        }
        break;
    }
    return V;
}

//-----------------------------------------------------------------------------
//! This function calculates the (initial) volume of an element. In some case, the volume
//! may only be approximate.
double FEMesh::CurrentElementVolume(RgElement& el)
{
    double V = 0;
    switch (el.elementType())
    {
        case FE_ELEM_SOLID:
        {
            RgSolidDomain* dom = dynamic_cast<RgSolidDomain*>(el.getDomain());
            assert(dom);
            // if (dom) V =dom->CurrentVolume(static_cast<RgSolidElement&>(el));
        }
        break;
        case FE_ELEM_SHELL:
        {
            /*FEShellDomain* dom = dynamic_cast<FEShellDomain*>(el.GetMeshPartition()); assert(dom);
            if (dom) return dom->CurrentVolume(static_cast<FEShellElement&>(el));*/
        }
        break;
    }
    return V;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by name

FENodeSet* FEMesh::FindNodeSet(const std::string& name)
{
    for (size_t i = 0; i < m_NodeSet.size(); ++i)
        if (m_NodeSet[i]->GetName() == name)
            return m_NodeSet[i];
    return 0;
}

//-----------------------------------------------------------------------------
//! Find a segment set set by name

// FESegmentSet* FEMesh::FindSegmentSet(const std::string& name)
//{
//	for (size_t i=0; i<m_LineSet.size(); ++i) if (m_LineSet[i]->GetName() ==  name) return m_LineSet[i];
//	return 0;
// }

//-----------------------------------------------------------------------------
//! Find a surface set set by name

FESurface* FEMesh::FindSurface(const std::string& name)
{
    /*for (size_t i = 0; i < m_Surf.size(); ++i)
        if (m_Surf[i]->GetName() == name)
            return m_Surf[i];*/
    return 0;
}

//-----------------------------------------------------------------------------
//! Find a surface set set by name and returns its index

int FEMesh::FindSurfaceIndex(const std::string& name)
{
    /*for (size_t i = 0; i < m_Surf.size(); ++i)
        if (m_Surf[i]->GetName() == name)
            return i;*/
    return -1;
}

FESurface* FEMesh::CreateSurface(FEFacetSet& facetSet)
{
    /* FESurface* surf = RANGO_NEW<FESurface>(m_fem, "");
     surf->Create(facetSet);
     AddSurface(surf);*/
    return nullptr;
}

//-----------------------------------------------------------------------------
//! Find a discrete element set set by name

// FEDiscreteSet* FEMesh::FindDiscreteSet(const std::string& name)
//{
//	for (size_t i=0; i<m_DiscSet.size(); ++i) if (m_DiscSet[i]->GetName() == name) return m_DiscSet[i];
//	return 0;
// }

//-----------------------------------------------------------------------------
//! Find a element set by name

RgElementSet* FEMesh::FindElementSet(const std::string& name)
{
    for (size_t i = 0; i < m_ElemSet.size(); ++i)
        if (m_ElemSet[i]->GetName() == name)
            return m_ElemSet[i];
    return 0;
}

//-----------------------------------------------------------------------------
FEFacetSet* FEMesh::FindFacetSet(const std::string& name)
{
    for (size_t i = 0; i < (int)m_FaceSet.size(); ++i)
        if (m_FaceSet[i]->GetName() == name)
            return m_FaceSet[i];
    return 0;
}

//-----------------------------------------------------------------------------
FESurfacePair* FEMesh::FindSurfacePair(const std::string& name)
{
    for (size_t i = 0; i < m_SurfPair.size(); ++i)
        if (m_SurfPair[i]->GetName() == name)
            return m_SurfPair[i];
    return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::Domains() const
{
    return (int)m_Domain.size();
}

//-----------------------------------------------------------------------------
RgDomain& FEMesh::Domain(int n)
{
    return *m_Domain[n];
}

FEModel* FEMesh::GetFEModel() const
{
    return m_fem;
}

void FEMesh::AddNodeSet(FENodeSet* pns)
{
    m_NodeSet.push_back(pns);
}

FENodeSet* FEMesh::NodeSet(int i)
{
    return m_NodeSet[i];
}

void FEMesh::AddElementSet(RgElementSet* pg)
{
    m_ElemSet.push_back(pg);
}

RgElementSet& FEMesh::ElementSet(int n)
{
    return *m_ElemSet[n];
}

FESurface& FEMesh::Surface(int n)
{
    return *m_Surf[n];
}

void FEMesh::AddSurface(FESurface* ps)
{
    m_Surf.push_back(ps);
}

void FEMesh::AddEdge(FEEdge* ps)
{
    m_Edge.push_back(ps);
}

FEFacetSet& FEMesh::FacetSet(int n)
{
    return *m_FaceSet[n];
}

void FEMesh::AddFacetSet(FEFacetSet* ps)
{
    m_FaceSet.push_back(ps);
}

FESurfacePair& FEMesh::SurfacePair(int n)
{
    return *m_SurfPair[n];
}

void FEMesh::AddSurfacePair(FESurfacePair* ps)
{
    m_SurfPair.push_back(ps);
}

FENodeElemList& FEMesh::NodeElementList()
{
    if (m_NEL.Size() != m_Node.size())
    {
        //m_NEL.Create(*this);  why link error
    }
        
    return m_NEL;
}

//-----------------------------------------------------------------------------
void FEMesh::AddDomain(RgDomain* pd)
{
    int N = (int)m_Domain.size();
    pd->SetID(N);
    m_Domain.push_back(pd);
    if (m_LUT)
        delete m_LUT;
    m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Find a domain

RgDomain* FEMesh::FindDomain(const std::string& name)
{
    for (size_t i = 0; i < m_Domain.size(); ++i)
        if (m_Domain[i]->GetName() == name)
            return m_Domain[i];
    return 0;
}

int FEMesh::FindDomainIndex(const std::string& name)
{
    for (size_t i = 0; i < m_Domain.size(); ++i)
        if (m_Domain[i]->GetName() == name)
            return i;
    return -1;
}

RgDomain* FEMesh::FindDomain(int domId)
{
    for (size_t i = 0; i < m_Domain.size(); ++i)
        if (m_Domain[i]->GetID() == domId)
            return m_Domain[i];
    return 0;
}

//-----------------------------------------------------------------------------
//! return an element
RgElement* FEMesh::Element(int n)
{
    if (n < 0)
        return nullptr;
    for (int i = 0; i < Domains(); ++i)
    {
        RgDomain& dom = Domain(i);
        int NEL = dom.Elements();
        if (n < NEL)
            return &dom.ElementRef(n);
        else
            n -= NEL;
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
//! Find a node from a given ID. return 0 if the node cannot be found.

FENode* FEMesh::FindNodeFromID(int nid)
{
    for (int i = 0; i < Nodes(); ++i)
    {
        FENode& node = Node(i);
        if (node.GetID() == nid)
            return &node;
    }
    return 0;
}

//-----------------------------------------------------------------------------
//! Find an element from a given ID. return 0 if the element cannot be found.

RgElement* FEMesh::FindElementFromID(int nid)
{
    if (m_LUT == 0)
        m_LUT = new FEElementLUT(*this);
    return m_LUT->Find(nid);
}

/*
RgElement* FEMesh::FindElementFromID(int nid)
{
    RgElement* pe = 0;

    for (int i=0; i<Domains(); ++i)
    {
        RgDomain& d = Domain(i);
        pe = d.FindElementFromID(nid);
        if (pe) return pe;
    }

    return pe;
}
*/


//-----------------------------------------------------------------------------
//! See if all elements are of a particular shape
bool FEMesh::IsType(ElementShape eshape)
{
    FEElementList elemList(*this);
    for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it)
    {
        RgElement& el = *it;
        if (el.elementShape() != eshape)
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
// Find the element in which point y lies
// RgSolidElement* FEMesh::FindSolidElement(Vector3d y, double r[3])
//{
//	int ND = (int) m_Domain.size();
//	for (int i=0; i<ND; ++i)
//	{
//		if (m_Domain[i]->Class() == FE_DOMAIN_SOLID)
//		{
//			RgSolidDomain& bd = static_cast<RgSolidDomain&>(*m_Domain[i]);
//			RgSolidElement* pe = bd.FindElement(y, r);
//			if (pe) return pe;
//		}
//	}
//	return 0;
//}

//-----------------------------------------------------------------------------
void FEMesh::ClearDomains()
{
    int N = Domains();
    for (int i = 0; i < N; ++i)
        delete m_Domain[i];
    m_Domain.clear();
    if (m_LUT)
        delete m_LUT;
    m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Rebuild the LUT
void FEMesh::RebuildLUT()
{
    if (m_LUT)
        delete m_LUT;
    m_LUT = new FEElementLUT(*this);
}

//-----------------------------------------------------------------------------
//! Calculate the surface representing the element boundaries
//! boutside : include all exterior facets
//! binside  : include all interior facets
FESurface* FEMesh::ElementBoundarySurface(bool boutside, bool binside)
{
    return nullptr;
    // if ((boutside == false) && (binside == false)) return 0;

    //// create the element neighbor list
    // FEElemElemList EEL;
    // EEL.Create(this);

    //// get the number of elements in this mesh
    // int NE = Elements();

    //// count the number of facets we have to create
    // int NF = 0;
    // FEElementList EL(*this);
    // FEElementList::iterator it = EL.begin();
    // for (int i=0; i<NE; ++i, ++it)
    //{
    //	RgElement& el = *it;
    //	int nf = el.Faces();
    //	for (int j=0; j<nf; ++j)
    //	{
    //		RgElement* pen = EEL.Neighbor(i, j);
    //		if ((pen == 0) && boutside) ++NF;
    //		if ((pen != 0) && (el.getId() < pen->getId()) && binside ) ++NF;
    //	}
    // }
    //// create the surface
    // FESurface* ps = RANGO_NEW<FESurface>( GetFEModel(),"");
    // if (NF == 0) return ps;
    // ps->Create(NF);

    //// build the surface elements
    // int face[RgElement::MAX_NODES];
    // NF = 0;
    // it = EL.begin();
    // for (int i=0; i<NE; ++i, ++it)
    //{
    //	RgElement& el = *it;
    //	int nf = el.Faces();
    //	for (int j=0; j<nf; ++j)
    //	{
    //		RgElement* pen = EEL.Neighbor(i, j);
    //		if (((pen == 0) && boutside)||
    //			((pen != 0) && (el.getId() < pen->getId()) && binside ))
    //		{
    //			FESurfaceElement& se = ps->Element(NF++);
    //			int faceNodes = el.GetFace(j, face);

    //			switch (faceNodes)
    //			{
    //			case 4: se.setType(FE_QUAD4G4); break;
    //			case 8: se.setType(FE_QUAD8G9); break;
    //			case 9: se.setType(FE_QUAD9G9); break;
    //			case 3: se.setType(FE_TRI3G1 ); break;
    //			case 6: se.setType(FE_TRI6G7 ); break;
    //			case 7: se.setType(FE_TRI7G7); break;
    //			default:
    //				assert(false);
    //			}
    //
    //			se.m_elem[0] = &el;
    //			if (pen) se.m_elem[1] = pen;
    //
    //			int nn = se.NodeSize();
    //			for (int k=0; k<nn; ++k)
    //			{
    //				//Rango TODO:
    //				//se.setNode(k) = face[k];
    //			}
    //		}
    //	}
    //}

    //// initialize the surface.
    //// This will set the local surface element ID's and also set the m_nelem IDs.
    // ps->Init();

    //// all done
    // return ps;
}

FESurface* FEMesh::ElementBoundarySurface(std::vector<RgDomain*> domains, bool boutside, bool binside)
{
    return nullptr;
    // if ((boutside == false) && (binside == false)) return nullptr;

    //// create the element neighbor list
    // FEElemElemList EEL;
    // EEL.Create(this);

    //// get the number of elements in this mesh
    // int NE = Elements();

    //// count the number of facets we have to create
    // int NF = 0;

    // for (int i = 0; i < domains.size(); i++)
    //{
    //	for (int j = 0; j < domains[i]->Elements(); j++)
    //	{
    //		RgElement& el = domains[i]->ElementRef(j);
    //		int nf = el.Faces();
    //		for (int k = 0; k<nf; ++k)
    //		{
    //			RgElement* pen = EEL.Neighbor(el.getId()-1, k);
    //			if ((pen == nullptr) && boutside) ++NF;
    //			else if (pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) &&
    //boutside) ++NF; 			if ((pen != nullptr) && (el.getId() < pen->getId()) && binside && (std::find(domains.begin(),
    //domains.end(), pen->GetMeshPartition()) != domains.end())) ++NF;
    //		}
    //	}
    // }

    //// create the surface
    // FESurface* ps = RANGO_NEW<FESurface>( GetFEModel(),"");
    // if (NF == 0) return ps;
    // ps->Create(NF);

    //// build the surface elements
    // int face[RgElement::MAX_NODES];
    // NF = 0;
    // for (int i = 0; i < domains.size(); i++)
    //{
    //	for (int j = 0; j < domains[i]->Elements(); j++)
    //	{
    //		RgElement& el = domains[i]->ElementRef(j);
    //		int nf = el.Faces();
    //		for (int k = 0; k < nf; ++k)
    //		{
    //			RgElement* pen = EEL.Neighbor(el.getId()-1, k);
    //			if (((pen == nullptr) && boutside) ||
    //				(pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) && boutside)
    //||
    //				((pen != nullptr) && (el.getId() < pen->getId()) && binside && (std::find(domains.begin(), domains.end(),
    //pen->GetMeshPartition()) != domains.end())))
    //			{
    //				FESurfaceElement& se = ps->Element(NF++);
    //				int faceNodes = el.GetFace(k, face);

    //				switch (faceNodes)
    //				{
    //				case 4: se.setType(FE_QUAD4G4); break;
    //				case 8: se.setType(FE_QUAD8G9); break;
    //				case 9: se.setType(FE_QUAD9G9); break;
    //				case 3: se.setType(FE_TRI3G1 ); break;
    //				case 6: se.setType(FE_TRI6G7 ); break;
    //				case 7: se.setType(FE_TRI7G7 ); break;
    //				default:
    //					assert(false);
    //				}

    //				se.m_elem[0] = &el;
    //				if (pen) se.m_elem[1] = pen;

    //				int nn = se.NodeSize();
    //				for (int p = 0; p < nn; ++p)
    //				{
    //					//Rango TODO:
    //					//se.m_node[p] = face[p];
    //				}
    //			}
    //		}
    //	}
    //}

    //// initialize the surface.
    //// This will set the local surface element ID's and also set the m_nelem IDs.
    // ps->Init();

    //// all done
    // return ps;
}

FEFacetSet* FEMesh::DomainBoundary(std::vector<RgDomain*> domains, bool boutside, bool binside)
{
    return nullptr;
    // if ((boutside == false) && (binside == false)) return nullptr;

    //// create the element neighbor list
    // FEElemElemList EEL;
    // EEL.Create(this);

    //// get the number of elements in this mesh
    // int NE = Elements();

    //// count the number of facets we have to create
    // int NF = 0;

    // for (int i = 0; i < domains.size(); i++)
    //{
    //	for (int j = 0; j < domains[i]->Elements(); j++)
    //	{
    //		RgElement& el = domains[i]->ElementRef(j);
    //		int nf = el.Faces();
    //		for (int k = 0; k < nf; ++k)
    //		{
    //			RgElement* pen = EEL.Neighbor(el.getId() - 1, k);
    //			if ((pen == nullptr) && boutside) ++NF;
    //			else if (pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) &&
    //boutside) ++NF; 			if ((pen != nullptr) && (el.getId() < pen->getId()) && binside && (std::find(domains.begin(),
    //domains.end(), pen->GetMeshPartition()) != domains.end())) ++NF;
    //		}
    //	}
    // }

    //// create the surface
    // FEFacetSet* ps = new FEFacetSet(GetFEModel());
    // if (NF == 0) return ps;
    // ps->Create(NF);

    //// build the surface elements
    // int faceNodes[RgElement::MAX_NODES];
    // NF = 0;
    // for (int i = 0; i < domains.size(); i++)
    //{
    //	for (int j = 0; j < domains[i]->Elements(); j++)
    //	{
    //		RgElement& el = domains[i]->ElementRef(j);
    //		int nf = el.Faces();
    //		for (int k = 0; k < nf; ++k)
    //		{
    //			RgElement* pen = EEL.Neighbor(el.getId() - 1, k);
    //			if (((pen == nullptr) && boutside) ||
    //				(pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) && boutside)
    //||
    //				((pen != nullptr) && (el.getId() < pen->getId()) && binside && (std::find(domains.begin(), domains.end(),
    //pen->GetMeshPartition()) != domains.end())))
    //			{
    //				FEFacetSet::FACET& f = ps->Face(NF++);
    //				int fn = el.GetFace(k, faceNodes);

    //				switch (fn)
    //				{
    //				case 4: f.ntype = FEFacetSet::FACET::QUAD4; break;
    //				case 8: f.ntype = FEFacetSet::FACET::QUAD8; break;
    //				case 9: f.ntype = FEFacetSet::FACET::QUAD9; break;
    //				case 3: f.ntype = FEFacetSet::FACET::TRI3; break;
    //				case 6: f.ntype = FEFacetSet::FACET::TRI6; break;
    //				case 7: f.ntype = FEFacetSet::FACET::TRI7; break;
    //				default:
    //					assert(false);
    //				}

    //				for (int p = 0; p < fn; ++p)
    //				{
    //					f.node[p] = faceNodes[p];
    //				}
    //			}
    //		}
    //	}
    //}

    //// all done
    // return ps;
}

//-----------------------------------------------------------------------------
//! Retrieve the nodal coordinates of an element in the reference configuration.
void FEMesh::GetInitialNodalCoordinates(const RgElement& el, Vector3d* node)
{
   /* const int neln = el.NodeSize();
    for (int i = 0; i < neln; ++i)
        node[i] = Node(el.m_node[i]).m_r0;*/
}

//-----------------------------------------------------------------------------
//! Retrieve the nodal coordinates of an element in the current configuration.
void FEMesh::GetNodalCoordinates(const RgElement& el, Vector3d* node)
{
   /* const int neln = el.NodeSize();
    for (int i = 0; i < neln; ++i)
        node[i] = Node(el.m_node[i]).m_rt;*/
}

//=============================================================================
FEElementLUT::FEElementLUT(FEMesh& mesh)
{
    // get the ID ranges
    m_minID = -1;
    m_maxID = -1;
    int NDOM = mesh.Domains();
    for (int i = 0; i < NDOM; ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        int NE = dom.Elements();
        for (int j = 0; j < NE; ++j)
        {
            RgElement& el = dom.ElementRef(j);
            int eid = el.getId();
            if ((eid < m_minID) || (m_minID == -1))
                m_minID = eid;
            if ((eid > m_maxID) || (m_maxID == -1))
                m_maxID = eid;
        }
    }

    // allocate size
    int nsize = m_maxID - m_minID + 1;
    m_elem.resize(nsize, (RgElement*)0);

    // fill the table
    for (int i = 0; i < NDOM; ++i)
    {
        RgDomain& dom = mesh.Domain(i);
        int NE = dom.Elements();
        for (int j = 0; j < NE; ++j)
        {
            RgElement& el = dom.ElementRef(j);
            int eid = el.getId();
            m_elem[eid - m_minID] = &el;
        }
    }
}

// Find an element from its ID
RgElement* FEElementLUT::Find(int nid)
{
    if ((nid < m_minID) || (nid > m_maxID))
        return 0;
    return m_elem[nid - m_minID];
}

// update the domains of the mesh
void FEMesh::Update(const FETimeInfo& tp)
{
    for (int i = 0; i < Domains(); ++i)
    {
        RgDomain& dom = Domain(i);
        // if (dom.IsActive()) dom.Update(tp);
    }
}
