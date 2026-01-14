/*********************************************************************
 * \file   FEMesh.h
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "datastructure/FEBoundingBox.h"
#include "FEDiscreteSet.h"
#include "elements/RgElementSet.h"
#include "FEFacetSet.h"
#include "FENode.h"
#include "FENodeElemList.h"
#include "femcore/FENodeSet.h"
//#include "elements/FEShellElement.h"
#include "elements/RgElementSet.h"
#include "FESurfacePair.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEEdge;
class FESurface;
class RgDomain;
class FEModel;
class FETimeInfo;
class FEDataMap;
class DumpStream;
class RgElement;

//---------------------------------------------------------------------------------------
// Helper class for faster lookup of elements based on their ID
class FEM_EXPORT FEElementLUT
{
public:
    FEElementLUT(FEMesh& mesh);

    // Find an element from its ID
    RgElement* Find(int nid);

private:
    std::vector<RgElement*> m_elem;
    int m_minID, m_maxID;
};

//-----------------------------------------------------------------------------
//! Defines a finite element mesh

//! All the geometry data is stored in this class.

class FEM_EXPORT FEMesh
{
public:
    //! constructor
    FEMesh(FEModel* fem);

    //! destructor
    virtual ~FEMesh();

    //! initialize material points in mesh
    void InitMaterialPoints();

    //! clear the mesh
    void Clear();

    //! allocate storage for mesh data
    void CreateNodes(int nodes);
    void AddNodes(int nodes);

    //! return number of nodes
    int Nodes() const;

    //! return total nr of elements
    int Elements() const;

    //! return the nr of elements of a specific domain type
    int Elements(int ndom_type) const;

    //! return reference to a node
    FENode& Node(int i);
    const FENode& Node(int i) const;

    const std::vector<FENode>& AllNode() const;

    //! Set the number of degrees of freedom on this mesh
    void SetDOFS(int n);

    //! update bounding box
    void UpdateBox();

    //! retrieve the bounding box
    FEBoundingBox& GetBoundingBox();

    //! remove isolated vertices
    int RemoveIsolatedVertices();

    //! Reset the mesh data
    void Reset();

    //! Calculates an elements volume in reference configuration
    double ElementVolume(RgElement& el);

    //! calculates element volume in current configuration
    double CurrentElementVolume(RgElement& el);

    //! Finds a node from a given ID
    FENode* FindNodeFromID(int nid);

    //! return an element (expensive way!)
    RgElement* Element(int i);

    //! Finds an element from a given ID
    RgElement* FindElementFromID(int nid);

    FENodeElemList& NodeElementList();

    //! See if all elements are of a particular shape
    bool IsType(ElementShape eshape);

    // --- NODESETS ---
    //! adds a node set to the mesh
    void AddNodeSet(FENodeSet* pns);

    //! number of nodesets
    int NodeSets() const;

    //! return a node set
    FENodeSet* NodeSet(int i);

    //! Find a nodeset by name
    FENodeSet* FindNodeSet(const std::string& name);

    // --- ELEMENT SETS ---
    int ElementSets() const;
    RgElementSet& ElementSet(int n);
    void AddElementSet(RgElementSet* pg);

    //! Find a element set by name
    RgElementSet* FindElementSet(const std::string& name);

    // --- DOMAINS ---
    int Domains() const;
    RgDomain& Domain(int n);
    const std::vector<RgDomain*>& AllDomain() const;

    void AddDomain(RgDomain* pd);

    RgDomain* FindDomain(const std::string& name);
    int FindDomainIndex(const std::string& name);
    RgDomain* FindDomain(int domId);

    //! clear all domains
    void ClearDomains();

    //! Rebuild the LUT
    void RebuildLUT();

    // --- SURFACES ---
    int Surfaces() const;
    FESurface& Surface(int n);
    void AddSurface(FESurface* ps);
    FESurface* FindSurface(const std::string& name);
    int FindSurfaceIndex(const std::string& name);

    // create a surface from a facet set
    FESurface* CreateSurface(FEFacetSet& facetSet);

    // --- EDGES ---
    int Edges() const;
    FEEdge& Edge(int n);
    void AddEdge(FEEdge* ps);

    // --- FACETSETS ---
    int FacetSets() const;
    FEFacetSet& FacetSet(int n);
    void AddFacetSet(FEFacetSet* ps);
    FEFacetSet* FindFacetSet(const std::string& name);

    // --- surface pairs ---
    int SurfacePairs() const;
    FESurfacePair& SurfacePair(int n);
    void AddSurfacePair(FESurfacePair* ps);
    FESurfacePair* FindSurfacePair(const std::string& name);

public:
    //! stream mesh data
    void Serialize(DumpStream& dmp);

    static void SaveClass(DumpStream& ar, FEMesh* p);
    static FEMesh* LoadClass(DumpStream& ar, FEMesh* p);

    // create a copy of this mesh
    void CopyFrom(FEMesh& mesh);

public:
    //! Calculate the surface representing the element boundaries
    //! boutside : include all exterior facets
    //! binside  : include all interior facets
    FESurface* ElementBoundarySurface(bool boutside = true, bool binside = false);

    //! Calculate the surface representing the element boundaries
    //! domains  : a list of which domains to create the surface from
    //! boutside : include all exterior facets
    //! binside  : include all interior facetsk
    FESurface* ElementBoundarySurface(std::vector<RgDomain*> domains, bool boutside = true, bool binside = false);
    FEFacetSet* DomainBoundary(std::vector<RgDomain*> domains, bool boutside = true, bool binside = false);

    //! get the nodal coordinates in reference configuration
    void GetInitialNodalCoordinates(const RgElement& el, Vector3d* node);

    //! get the nodal coordinates in current configuration
    void GetNodalCoordinates(const RgElement& el, Vector3d* node);

    // Get the FE model
    FEModel* GetFEModel() const;

    // update the domains of the mesh
    void Update(const FETimeInfo& tp);

public:  // data maps
    void ClearDataMaps();
    void AddDataMap(FEDataMap* map);
    FEDataMap* FindDataMap(const std::string& map);

    int DataMaps() const;
    FEDataMap* GetDataMap(int i);

private:
    std::vector<FENode> m_Node;              //!< nodes
    std::vector<RgDomain*> m_Domain;         //!< list of domains
    std::vector<FESurface*> m_Surf;          //!< surfaces
    std::vector<FEEdge*> m_Edge;             //!< Edges

    std::vector<FENodeSet*> m_NodeSet;       //!< node sets
    std::vector<RgElementSet*> m_ElemSet;    //!< element sets
};