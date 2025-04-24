#pragma once
#include "femcore/fem_export.h"

#include <string>

class FEMesh;
class DumpStream;

class FEM_EXPORT FEItemList
{
public:
    FEItemList(FEMesh* mesh);
    virtual ~FEItemList();

    // get the mesh
    FEMesh* GetMesh() const;

    void SetMesh(FEMesh* mesh);

    const std::string& GetName() const;
    void SetName(const std::string& name);

public:
    void Serialize(DumpStream& ar);

    static FEItemList* LoadClass(DumpStream& ar, FEItemList* p);
    static void SaveClass(DumpStream& ar, FEItemList* p);

protected:
    FEMesh* m_mesh;

    std::string m_name;
};
