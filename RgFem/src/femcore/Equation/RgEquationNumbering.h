/*********************************************************************
 * \file   FEEquationNumbering.h
 * \brief  Global equation numbering system for FE model
 *
 * This class handles the numbering of equations for the global system,
 * taking into account DOF activation and boundary conditions.
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#pragma once

#include "femcore/fem_export.h"

#include <map>
#include <vector>

class FEModel;
class FEMesh;
class FENode;
class RgDofSchema;

//=============================================================================
// FEEquationNumbering - 全局方程编号管理器
//=============================================================================

/// 管理全局方程编号
class FEM_EXPORT RgEquationNumbering
{
public:
    RgEquationNumbering(FEModel* model);
    ~RgEquationNumbering();

    //=========================================================================
    // 方程编号
    //=========================================================================

    /// 对所有节点DOF进行全局编号
    /// @return 总方程数
    int NumberEquations();

    /// 对指定节点集的DOF进行编号
    /// @param nodeIds 节点ID列表
    /// @return 为这些节点分配的方程数
    int NumberNodeSet(const std::vector<int>& nodeIds);

    /// 重新编号（当边界条件改变后）
    int RenumberEquations();

    //=========================================================================
    // 查询接口
    //=========================================================================

    int GetTotalEquations() const
    {
        return m_neq;
    }

    int GetFreeDofCount() const
    {
        return m_nFreeDofs;
    }

    int GetPrescribedDofCount() const
    {
        return m_nPrescribedDofs;
    }

    /// 获取未激活DOF数量
    int GetInactiveDofCount() const
    {
        return m_nInactiveDofs;
    }

    /// 根据节点ID和DOF索引获取方程号
    int GetEquation(int nodeId, int dofIndex) const;

    //=========================================================================
    // 统计和诊断
    //=========================================================================

    /// 打印编号统计信息
    void PrintStatistics() const;

    /// 打印详细的DOF分布
    void PrintDofDistribution() const;

    /// 验证编号的正确性
    bool Validate() const;

    /// 导出编号信息到文件
    bool ExportNumbering(const char* filename) const;

private:
    //=========================================================================
    // 内部方法
    //=========================================================================

    /// 重置所有方程号
    void ResetEquations();

    /// 收集统计信息
    void CollectStatistics();

    /// 按DOF类型统计
    void CollectDofTypeStatistics();

private:
    FEModel* m_model;                        ///< 指向FE模型

    int m_neq;                               ///< 总方程数
    int m_nFreeDofs;                         ///< 自由DOF数量
    int m_nPrescribedDofs;                   ///< 固定DOF数量
    int m_nInactiveDofs;                     ///< 未激活DOF数量

    std::map<int, int> m_dofTypeFree;        ///< 每种DOF的自由数量
    std::map<int, int> m_dofTypePrescribed;  ///< 每种DOF的固定数量
    std::map<int, int> m_dofTypeInactive;    ///< 每种DOF的未激活数量
};

//=============================================================================
// 辅助函数
//=============================================================================

namespace FEEquationUtil
{
    /// 从节点和DOF索引计算紧凑的全局DOF编号（仅用于调试）
    FEM_EXPORT int CompactDofId(int nodeId, int dofIndex, int dofsPerNode);

    /// 从紧凑编号解析出节点ID和DOF索引
    FEM_EXPORT void ParseCompactDofId(int compactId, int dofsPerNode, int& nodeId, int& dofIndex);
}  // namespace FEEquationUtil
