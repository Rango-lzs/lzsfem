#pragma once
#include "femcore/FEItemList.h"
#include "femcore/fem_export.h"
#include "femcore/FENodeList.h"

#include <string>
#include <vector>

//-----------------------------------------------------------------------------
// Forward declarations
class FEMesh;
class FENode;
class DumpStream;
class FEModel;

//-----------------------------------------------------------------------------
//! Defines a node set of the model
//
class FEM_EXPORT FENodeSet : public FEItemList
{
public:
    FENodeSet(FEModel* fem);

    void Add(int n);
    void Add(const std::vector<int>& ns);
    void Add(const FENodeList& nodeList);

    void Clear();

    int Size() const
    {
        return m_Node.Size();
    }

    int operator[](int i) const
    {
        return m_Node[i];
    }

    FENodeList GetNodeList()
    {
        return m_Node;
    }
    const FENodeList& GetNodeList() const
    {
        return m_Node;
    }

    FENode* Node(int i);
    const FENode* Node(int i) const;

public:
    void Serialize(DumpStream& ar);
    static void SaveClass(DumpStream& ar, FENodeSet* p);
    static FENodeSet* LoadClass(DumpStream& ar, FENodeSet* p);

protected:
    FENodeList m_Node;  //!< list of nodes
};
