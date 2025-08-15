// RgElement.h
#pragma once
#include "ElementShape/ElementShape.h"  // 假设存在形状函数定义
#include "materials/Material.h"        // 材料基类
#include "datastructure/Matrix.h"      // 矩阵工具

class RgElement {
public:
    // 构造/析构函数
    RgElement(int id, ElementShape* shape, Material* mat);
    virtual ~RgElement();

    // 核心接口
    virtual void ComputeStiffnessMatrix(Matrix& ke) = 0;  // 计算刚度矩阵
    virtual void ComputeInternalForces(Vector& fe) = 0;   // 计算内力向量
    virtual void UpdateState() = 0;                       // 更新单元状态

    // 辅助接口
    int GetID() const { return m_id; }
    ElementShape* GetShape() const { return m_shape; }
    Material* GetMaterial() const { return m_mat; }
    
    // 积分点操作
    virtual int GetNumIntegrationPoints() const = 0;
    virtual const GaussPoint& GetIntegrationPoint(int n) const = 0;

protected:
    int m_id;                // 单元ID
    ElementShape* m_shape;   // 单元形状定义
    Material* m_mat;         // 材料属性
};