#pragma once
#include "femcore/fem_export.h"
#include <vector>

class RgMaterialPoint;

//记录单元材料点的状态数据
class FEM_EXPORT RgElementState
{
public:
    RgElementState();
    ~RgElementState();

    RgElementState(const RgElementState& s);
    RgElementState& operator=(const RgElementState& s);

    void Create(int n);
    RgMaterialPoint*& operator[](int i);
    const std::vector<RgMaterialPoint*>& getMatPoints() const;

private:
    void Clear();
    std::vector<RgMaterialPoint*> m_data;
};