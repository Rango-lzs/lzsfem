#pragma once
#include "femcore/fem_export.h"
#include <vector>

class FEMaterialPoint;

class FEM_EXPORT FEElementState
{
public:
    FEElementState();
    ~FEElementState();

    FEElementState(const FEElementState& s);
    FEElementState& operator=(const FEElementState& s);

    void Create(int n);
    FEMaterialPoint*& operator[](int i);
    const std::vector<FEMaterialPoint*>& getMatPoints() const;

private:
    void Clear();
    std::vector<FEMaterialPoint*> m_data;
};