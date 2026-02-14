/*********************************************************************
 * \file   FEEquationNumbering.cpp
 * \brief  Implementation of global equation numbering
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#include "femcore/Equation/RgEquationNumbering.h"

#include "femcore/FEMesh.h"
#include "femcore/FEModel.h"
#include "femcore/FENode.h"
#include "femcore/RgDofSchema.h"
#include "logger/log.h"

#include <algorithm>
#include <fstream>
#include <iomanip>

 //=============================================================================
 // FEEquationNumbering Implementation
 //=============================================================================

RgEquationNumbering::RgEquationNumbering(FEModel* model)
    : m_model(model)
    , m_neq(0)
    , m_nFreeDofs(0)
    , m_nPrescribedDofs(0)
    , m_nInactiveDofs(0)
{
}

RgEquationNumbering::~RgEquationNumbering()
{
}

//-----------------------------------------------------------------------------
// 主要编号函数
//-----------------------------------------------------------------------------

int RgEquationNumbering::NumberEquations()
{
    RgLog("\n");
    RgLog("=======================================================\n");
    RgLog("Global Equation Numbering\n");
    RgLog("=======================================================\n");

    // 获取网格和DOF schema
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    RgLog("Nodes:          %d\n", nNodes);
    RgLog("DOFs per node:  %d\n", dofsPerNode);
    RgLog("Total DOF slots: %d\n", nNodes * dofsPerNode);
    RgLog("\n");

    // 重置所有方程号
    ResetEquations();

    // 方程号计数器
    int eqn = 0;

    // 遍历所有节点
    for (int i = 0; i < nNodes; ++i)
    {
        FENode& node = mesh.Node(i);

        // 遍历节点的每个DOF
        for (int d = 0; d < dofsPerNode; ++d)
        {
            // 只为激活且自由的DOF分配方程号
            if (node.is_active(d) && node.IsDofFree(d))
            {
                node.setDofIdx(d, eqn);
                eqn++;
            }
            else
            {
                // 未激活或固定的DOF标记为-1（无方程）
                node.setDofIdx(d, -1);
            }
        }
    }

    // 保存总方程数
    m_neq = eqn;

    RgLog("Numbering complete.\n");
    RgLog("Total equations: %d\n\n", m_neq);

    // 收集统计信息
    CollectStatistics();

    // 打印统计
    PrintStatistics();

    // 验证编号
    if (!Validate())
    {
        RgLogError("Equation numbering validation FAILED!");
        return -1;
    }

    RgLog("=======================================================\n\n");

    return m_neq;
}

//-----------------------------------------------------------------------------
int RgEquationNumbering::NumberNodeSet(const std::vector<int>& nodeIds)
{
    const RgDofSchema& schema = m_model->GetDofSchema();
    FEMesh& mesh = m_model->GetMesh();

    int dofsPerNode = schema.GetDofsPerNode();
    int count = 0;

    for (int nodeId : nodeIds)
    {
        if (nodeId < 0 || nodeId >= mesh.Nodes())
            continue;

        FENode& node = mesh.Node(nodeId);

        for (int d = 0; d < dofsPerNode; ++d)
        {
            if (node.is_active(d) && node.IsDofFree(d))
            {
                count++;
            }
        }
    }

    return count;
}

//-----------------------------------------------------------------------------
int RgEquationNumbering::RenumberEquations()
{
    RgLog("Re-numbering equations...\n");
    return NumberEquations();
}

//-----------------------------------------------------------------------------
int RgEquationNumbering::GetEquation(int nodeId, int dofIndex) const
{
    FEMesh& mesh = m_model->GetMesh();

    if (nodeId < 0 || nodeId >= mesh.Nodes())
        return -1;

    const FENode& node = mesh.Node(nodeId);
    return node.getDofIdx(dofIndex);
}

//-----------------------------------------------------------------------------
// 内部方法
//-----------------------------------------------------------------------------

void RgEquationNumbering::ResetEquations()
{
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    // 清除所有节点的方程号
    for (int i = 0; i < nNodes; ++i)
    {
        FENode& node = mesh.Node(i);
        for (int d = 0; d < dofsPerNode; ++d)
        {
            node.setDofIdx(d, -1);
        }
    }

    // 重置统计
    m_neq = 0;
    m_nFreeDofs = 0;
    m_nPrescribedDofs = 0;
    m_nInactiveDofs = 0;
    m_dofTypeFree.clear();
    m_dofTypePrescribed.clear();
    m_dofTypeInactive.clear();
}

//-----------------------------------------------------------------------------
void RgEquationNumbering::CollectStatistics()
{
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    m_nFreeDofs = 0;
    m_nPrescribedDofs = 0;
    m_nInactiveDofs = 0;

    // 遍历所有节点统计
    for (int i = 0; i < nNodes; ++i)
    {
        const FENode& node = mesh.Node(i);

        for (int d = 0; d < dofsPerNode; ++d)
        {
            if (!node.is_active(d))
            {
                m_nInactiveDofs++;
            }
            else if (node.IsDofPrescribed(d))
            {
                m_nPrescribedDofs++;
            }
            else if (node.IsDofFree(d))
            {
                m_nFreeDofs++;
            }
        }
    }

    // 按DOF类型统计
    CollectDofTypeStatistics();
}

//-----------------------------------------------------------------------------
void RgEquationNumbering::CollectDofTypeStatistics()
{
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    m_dofTypeFree.clear();
    m_dofTypePrescribed.clear();
    m_dofTypeInactive.clear();

    for (int d = 0; d < dofsPerNode; ++d)
    {
        m_dofTypeFree[d] = 0;
        m_dofTypePrescribed[d] = 0;
        m_dofTypeInactive[d] = 0;
    }

    for (int i = 0; i < nNodes; ++i)
    {
        const FENode& node = mesh.Node(i);

        for (int d = 0; d < dofsPerNode; ++d)
        {
            if (!node.is_active(d))
            {
                m_dofTypeInactive[d]++;
            }
            else if (node.IsDofPrescribed(d))
            {
                m_dofTypePrescribed[d]++;
            }
            else if (node.IsDofFree(d))
            {
                m_dofTypeFree[d]++;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// 统计和诊断
//-----------------------------------------------------------------------------

void RgEquationNumbering::PrintStatistics() const
{
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();
    int totalDofs = nNodes * dofsPerNode;

    RgLog("-------------------------------------------------------\n");
    RgLog("Equation Numbering Statistics\n");
    RgLog("-------------------------------------------------------\n");

    RgLog("Total DOF slots:       %8d\n", totalDofs);
    RgLog("\n");

    int activeDofs = m_nFreeDofs + m_nPrescribedDofs;
    RgLog("Active DOFs:           %8d (%.1f%%)\n", activeDofs, 100.0 * activeDofs / totalDofs);
    RgLog("  Free (equations):    %8d (%.1f%%)\n", m_nFreeDofs, 100.0 * m_nFreeDofs / totalDofs);
    RgLog("  Prescribed (BC):     %8d (%.1f%%)\n", m_nPrescribedDofs, 100.0 * m_nPrescribedDofs / totalDofs);
    RgLog("\n");

    RgLog("Inactive DOFs:         %8d (%.1f%%)\n", m_nInactiveDofs, 100.0 * m_nInactiveDofs / totalDofs);
    RgLog("\n");

    RgLog("Total Equations:       %8d\n", m_neq);
    RgLog("-------------------------------------------------------\n");

    // 打印DOF分布
    PrintDofDistribution();
}

//-----------------------------------------------------------------------------
void RgEquationNumbering::PrintDofDistribution() const
{
    const RgDofSchema& schema = m_model->GetDofSchema();
    int dofsPerNode = schema.GetDofsPerNode();

    RgLog("\nDOF Distribution by Type:\n");
    RgLog("-------------------------------------------------------\n");
    RgLog("%-8s %10s %10s %10s %10s\n", "DOF", "Active", "Free", "Fixed", "Inactive");
    RgLog("-------------------------------------------------------\n");

    for (int d = 0; d < dofsPerNode; ++d)
    {
        DofType dofType = schema.GetDofType(d);
        std::string dofName = schema.GetDofName(dofType);

        int nFree = m_dofTypeFree.count(d) ? m_dofTypeFree.at(d) : 0;
        int nPrescribed = m_dofTypePrescribed.count(d) ? m_dofTypePrescribed.at(d) : 0;
        int nInactive = m_dofTypeInactive.count(d) ? m_dofTypeInactive.at(d) : 0;
        int nActive = nFree + nPrescribed;

        RgLog("%-8s %10d %10d %10d %10d\n", dofName.c_str(), nActive, nFree, nPrescribed, nInactive);
    }

    RgLog("-------------------------------------------------------\n\n");
}

//-----------------------------------------------------------------------------
bool RgEquationNumbering::Validate() const
{
    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    RgLog("Validating equation numbering...\n");

    bool valid = true;
    int errorCount = 0;

    // 检查方程号的使用情况
    std::vector<bool> eqnUsed(m_neq, false);
    std::vector<int> eqnNode(m_neq, -1);
    std::vector<int> eqnDof(m_neq, -1);

    for (int i = 0; i < nNodes; ++i)
    {
        const FENode& node = mesh.Node(i);

        for (int d = 0; d < dofsPerNode; ++d)
        {
            int eqn = node.getDofIdx(d);

            if (eqn >= 0)
            {
                // 检查方程号范围
                if (eqn >= m_neq)
                {
                    RgLogError("  Node %d DOF %d: equation %d out of range [0, %d)", i, d, eqn, m_neq);
                    valid = false;
                    errorCount++;
                    continue;
                }

                // 检查重复使用
                if (eqnUsed[eqn])
                {
                    RgLogError("  Equation %d used by both Node %d DOF %d and Node %d DOF %d", eqn, eqnNode[eqn],
                        eqnDof[eqn], i, d);
                    valid = false;
                    errorCount++;
                }
                else
                {
                    eqnUsed[eqn] = true;
                    eqnNode[eqn] = i;
                    eqnDof[eqn] = d;
                }

                // 检查DOF必须是激活且自由的
                if (!node.is_active(d))
                {
                    RgLogError("  Node %d DOF %d: has equation but not active", i, d);
                    valid = false;
                    errorCount++;
                }

                if (node.IsDofPrescribed(d))
                {
                    RgLogError("  Node %d DOF %d: has equation but is prescribed", i, d);
                    valid = false;
                    errorCount++;
                }
            }
            else  // eqn < 0 (无方程号)
            {
                // 自由且激活的DOF必须有方程号
                if (node.is_active(d) && node.IsDofFree(d))
                {
                    RgLogError("  Node %d DOF %d: is free but has no equation", i, d);
                    valid = false;
                    errorCount++;
                }
            }
        }
    }

    // 检查所有方程号都被使用
    for (int e = 0; e < m_neq; ++e)
    {
        if (!eqnUsed[e])
        {
            RgLogError("  Equation %d is not assigned", e);
            valid = false;
            errorCount++;
        }
    }

    // 检查统计一致性
    if (m_nFreeDofs != m_neq)
    {
        RgLogError("  Inconsistency: %d free DOFs but %d equations", m_nFreeDofs, m_neq);
        valid = false;
        errorCount++;
    }

    // 报告结果
    if (valid)
    {
        RgLog("Validation PASSED\n\n");
    }
    else
    {
        RgLog("Validation FAILED with %d errors\n\n", errorCount);
    }

    return valid;
}

//-----------------------------------------------------------------------------
bool RgEquationNumbering::ExportNumbering(const char* filename) const
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        RgLogError("Cannot open file: %s", filename);
        return false;
    }

    FEMesh& mesh = m_model->GetMesh();
    const RgDofSchema& schema = m_model->GetDofSchema();

    int nNodes = mesh.Nodes();
    int dofsPerNode = schema.GetDofsPerNode();

    // 写入头部
    file << "======================================================\n";
    file << "Global Equation Numbering Report\n";
    file << "======================================================\n\n";

    // 统计信息
    file << "Summary:\n";
    file << "  Nodes:            " << nNodes << "\n";
    file << "  DOFs per node:    " << dofsPerNode << "\n";
    file << "  Total equations:  " << m_neq << "\n";
    file << "  Free DOFs:        " << m_nFreeDofs << "\n";
    file << "  Prescribed DOFs:  " << m_nPrescribedDofs << "\n";
    file << "  Inactive DOFs:    " << m_nInactiveDofs << "\n\n";

    // DOF类型分布
    file << "DOF Type Distribution:\n";
    file << std::left << std::setw(10) << "DOF" << std::setw(12) << "Active" << std::setw(12) << "Free" << std::setw(12)
        << "Prescribed" << std::setw(12) << "Inactive"
        << "\n";
    file << std::string(58, '-') << "\n";

    for (int d = 0; d < dofsPerNode; ++d)
    {
        DofType dofType = schema.GetDofType(d);
        std::string dofName = schema.GetDofName(dofType);

        int nFree = m_dofTypeFree.count(d) ? m_dofTypeFree.at(d) : 0;
        int nPrescribed = m_dofTypePrescribed.count(d) ? m_dofTypePrescribed.at(d) : 0;
        int nInactive = m_dofTypeInactive.count(d) ? m_dofTypeInactive.at(d) : 0;
        int nActive = nFree + nPrescribed;

        file << std::left << std::setw(10) << dofName << std::setw(12) << nActive << std::setw(12) << nFree
            << std::setw(12) << nPrescribed << std::setw(12) << nInactive << "\n";
    }

    file << "\n";

    // 详细节点信息
    file << "Detailed Node Equations:\n";
    file << "======================================================\n\n";

    file << std::left << std::setw(8) << "Node" << std::setw(8) << "DOF" << std::setw(12) << "Active" << std::setw(14)
        << "State" << std::setw(10) << "Equation" << std::setw(14) << "Value"
        << "\n";
    file << std::string(66, '-') << "\n";

    for (int i = 0; i < nNodes; ++i)
    {
        const FENode& node = mesh.Node(i);

        for (int d = 0; d < dofsPerNode; ++d)
        {
            DofType dofType = schema.GetDofType(d);
            std::string dofName = schema.GetDofName(dofType);

            bool active = node.is_active(d);
            FENode::DofState state = node.getDofState(d);
            int eqn = node.getDofIdx(d);
            double value = node.get(d);

            const char* stateStr = "INACTIVE";
            if (state == DOF_OPEN)
                stateStr = "OPEN";
            else if (state == DOF_PRESCRIBED)
                stateStr = "PRESCRIBED";

            file << std::left << std::setw(8) << i << std::setw(8) << dofName << std::setw(12)
                << (active ? "Yes" : "No") << std::setw(14) << stateStr << std::setw(10) << eqn << std::setw(14)
                << std::scientific << std::setprecision(6) << value << "\n";
        }
    }

    file << "\n======================================================\n";

    file.close();

    RgLog("Numbering exported to: %s\n", filename);

    return true;
}

//=============================================================================
// 辅助函数实现
//=============================================================================

namespace FEEquationUtil
{
    int CompactDofId(int nodeId, int dofIndex, int dofsPerNode)
    {
        return nodeId * dofsPerNode + dofIndex;
    }

    void ParseCompactDofId(int compactId, int dofsPerNode, int& nodeId, int& dofIndex)
    {
        nodeId = compactId / dofsPerNode;
        dofIndex = compactId % dofsPerNode;
    }
}  // namespace FEEquationUtil
