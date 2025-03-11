#pragma once
#include "femcore/fem_export.h"

#include <vector>

//-----------------------------------------------------------------------------
class FEMesh;
class FENode;
class DumpStream;

//-----------------------------------------------------------------------------
// Defines an array of node indices.
class FEM_EXPORT FENodeList
{
public:
    FENodeList(FEMesh* mesh = nullptr);
    FENodeList(const FENodeList& nodeList);
    FENodeList& operator=(const FENodeList& nodeList);

    int Size() const;

    void Add(int n);
    void Add(const std::vector<int>& nodeList);
    void Add(const FENodeList& nodeList);

    void Clear();

    int operator[](int n) const
    {
        return m_nodes[n];
    }
    FENode* Node(int i);
    const FENode* Node(int i) const;

    void Serialize(DumpStream& ar);

    FEMesh* GetMesh()
    {
        return m_mesh;
    }

    // returns the local index from a global (mesh based) node index,
    // or -1 if the node is not part of the set.
    int GlobalToLocalID(int globalId) const;

private:
    FEMesh* m_mesh;
    std::vector<int> m_nodes;
};
