/*********************************************************************
 * \file   FEEquationNumbering_Examples.cpp
 * \brief  Usage examples for the equation numbering system
 *
 * \date   February 2026
 *********************************************************************/

#include "FEEquationNumbering.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENode.h"
#include "femcore/Domain/RgSolidDomain.h"
#include "RgDofSchema.h"

//=============================================================================
// 示例1: 基本使用 - 简单的Solid模型
//=============================================================================

void Example1_BasicSolidModel()
{
    RgLog("\n");
    RgLog("=========================================================\n");
    RgLog("Example 1: Basic Solid Model Equation Numbering\n");
    RgLog("=========================================================\n\n");
    
    // 创建模型
    FEModel model;
    FEMesh& mesh = model.GetMesh();
    
    // 创建一个简单的网格：8个节点，1个hex8单元
    mesh.CreateNodes(8);
    
    // 设置节点坐标
    mesh.Node(0).m_r0 = Vector3d(0, 0, 0);
    mesh.Node(1).m_r0 = Vector3d(1, 0, 0);
    mesh.Node(2).m_r0 = Vector3d(1, 1, 0);
    mesh.Node(3).m_r0 = Vector3d(0, 1, 0);
    mesh.Node(4).m_r0 = Vector3d(0, 0, 1);
    mesh.Node(5).m_r0 = Vector3d(1, 0, 1);
    mesh.Node(6).m_r0 = Vector3d(1, 1, 1);
    mesh.Node(7).m_r0 = Vector3d(0, 1, 1);
    
    // 创建域和单元
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(1, ElementType::FE_HEX8G8);
    dom->SetName("SolidDomain");
    
    // 设置单元节点连接性
    RgElement& elem = dom->ElementRef(0);
    for (int i = 0; i < 8; ++i)
    {
        elem.setNodeId(i, i);
    }
    
    mesh.AddDomain(dom);
    
    // 初始化模型（会自动配置DOF schema并激活节点DOF）
    model.Init();
    
    // 应用边界条件：固定节点0的所有位移
    RgDofSchema& schema = model.GetDofSchema();
    FENode& node0 = mesh.Node(0);
    
    int uIdx = schema.GetDofIndex("u");
    int vIdx = schema.GetDofIndex("v");
    int wIdx = schema.GetDofIndex("w");
    
    node0.SetDofState(uIdx, DOF_PRESCRIBED);
    node0.set(uIdx, 0.0);
    
    node0.SetDofState(vIdx, DOF_PRESCRIBED);
    node0.set(vIdx, 0.0);
    
    node0.SetDofState(wIdx, DOF_PRESCRIBED);
    node0.set(wIdx, 0.0);
    
    // 创建方程编号器
    FEEquationNumbering numbering(&model);
    
    // 执行方程编号
    int neq = numbering.NumberEquations();
    
    RgLog("Expected: 7 nodes * 3 DOFs = 21 equations\n");
    RgLog("Got:      %d equations\n\n", neq);
    
    // 验证
    numbering.Validate();
    
    // 导出详细信息
    numbering.ExportNumbering("example1_numbering.txt");
}

//=============================================================================
// 示例2: 混合单元模型 - Solid + Shell
//=============================================================================

void Example2_MixedSolidShellModel()
{
    RgLog("\n");
    RgLog("=========================================================\n");
    RgLog("Example 2: Mixed Solid+Shell Model\n");
    RgLog("=========================================================\n\n");
    
    FEModel model;
    FEMesh& mesh = model.GetMesh();
    
    // 创建节点（假设10个节点）
    mesh.CreateNodes(10);
    
    // ... 设置节点坐标 ...
    
    // 创建solid域（使用节点0-7）
    RgSolidDomain* solidDom = new RgSolidDomain(&model);
    solidDom->Create(1, ElementType::FE_HEX8G8);
    solidDom->SetName("SolidRegion");
    mesh.AddDomain(solidDom);
    
    // 创建shell域（使用节点8-9，假设有shell domain）
    // RgShellDomain* shellDom = new RgShellDomain(&model);
    // ... 设置shell单元 ...
    // mesh.AddDomain(shellDom);
    
    // 初始化模型
    model.Init();
    
    /* 
    此时的DOF配置：
    - 如果只有solid: 每个节点3个DOF (u,v,w)
    - 如果有shell:   每个节点6个DOF (u,v,w,Rx,Ry,Rz)
    
    节点激活情况：
    - 只被solid使用的节点: u,v,w激活，Rx,Ry,Rz未激活
    - 只被shell使用的节点: 全部6个DOF激活
    - 共享节点: 全部6个DOF激活
    */
    
    // 应用边界条件
    RgDofSchema& schema = model.GetDofSchema();
    
    // 固定节点0的位移
    for (int d = 0; d < 3; ++d)
    {
        mesh.Node(0).SetDofState(d, DOF_PRESCRIBED);
        mesh.Node(0).set(d, 0.0);
    }
    
    // 方程编号
    FEEquationNumbering numbering(&model);
    int neq = numbering.NumberEquations();
    
    RgLog("Total equations: %d\n\n", neq);
    
    // 打印详细信息
    numbering.PrintStatistics();
    numbering.PrintDofDistribution();
}

//=============================================================================
// 示例3: 复杂边界条件
//=============================================================================

void Example3_ComplexBoundaryConditions()
{
    RgLog("\n");
    RgLog("=========================================================\n");
    RgLog("Example 3: Complex Boundary Conditions\n");
    RgLog("=========================================================\n\n");
    
    FEModel model;
    FEMesh& mesh = model.GetMesh();
    
    // 创建网格
    mesh.CreateNodes(20);
    
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(2, ElementType::FE_HEX8G8);
    mesh.AddDomain(dom);
    
    // 初始化
    model.Init();
    
    RgDofSchema& schema = model.GetDofSchema();
    
    // 场景：悬臂梁
    // - 节点0-3: 固定端，固定所有位移
    // - 节点4-7: 对称边界，固定v=0
    // - 其他节点: 自由
    
    int uIdx = schema.GetDofIndex("u");
    int vIdx = schema.GetDofIndex("v");
    int wIdx = schema.GetDofIndex("w");
    
    // 固定端（节点0-3）
    for (int i = 0; i < 4; ++i)
    {
        FENode& node = mesh.Node(i);
        node.SetDofState(uIdx, DOF_PRESCRIBED);
        node.SetDofState(vIdx, DOF_PRESCRIBED);
        node.SetDofState(wIdx, DOF_PRESCRIBED);
        node.set(uIdx, 0.0);
        node.set(vIdx, 0.0);
        node.set(wIdx, 0.0);
    }
    
    // 对称边界（节点4-7）
    for (int i = 4; i < 8; ++i)
    {
        FENode& node = mesh.Node(i);
        node.SetDofState(vIdx, DOF_PRESCRIBED);
        node.set(vIdx, 0.0);
    }
    
    // 方程编号
    FEEquationNumbering numbering(&model);
    int neq = numbering.NumberEquations();
    
    /*
    预期结果：
    - 节点0-3: 3个DOF固定 = 0个方程
    - 节点4-7: 1个DOF固定，2个自由 = 2个方程/节点
    - 节点8-19: 3个自由DOF = 3个方程/节点
    
    总方程数 = 0*4 + 2*4 + 3*12 = 8 + 36 = 44
    */
    
    RgLog("Expected equations: 44\n");
    RgLog("Got equations:      %d\n\n", neq);
    
    numbering.PrintStatistics();
    numbering.Validate();
}

//=============================================================================
// 示例4: 重新编号
//=============================================================================

void Example4_Renumbering()
{
    RgLog("\n");
    RgLog("=========================================================\n");
    RgLog("Example 4: Renumbering After Boundary Condition Change\n");
    RgLog("=========================================================\n\n");
    
    FEModel model;
    FEMesh& mesh = model.GetMesh();
    
    mesh.CreateNodes(10);
    
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(1, ElementType::FE_HEX8G8);
    mesh.AddDomain(dom);
    
    model.Init();
    
    RgDofSchema& schema = model.GetDofSchema();
    FEEquationNumbering numbering(&model);
    
    // 第一次编号：无边界条件
    RgLog("--- Initial numbering (no BCs) ---\n");
    int neq1 = numbering.NumberEquations();
    RgLog("Equations: %d\n\n", neq1);
    
    // 应用边界条件：固定节点0
    RgLog("--- Applying BC: Fix node 0 ---\n");
    int uIdx = schema.GetDofIndex("u");
    int vIdx = schema.GetDofIndex("v");
    int wIdx = schema.GetDofIndex("w");
    
    mesh.Node(0).SetDofState(uIdx, DOF_PRESCRIBED);
    mesh.Node(0).SetDofState(vIdx, DOF_PRESCRIBED);
    mesh.Node(0).SetDofState(wIdx, DOF_PRESCRIBED);
    
    // 重新编号
    RgLog("--- Renumbering ---\n");
    int neq2 = numbering.RenumberEquations();
    RgLog("Equations: %d\n\n", neq2);
    
    RgLog("Equation reduction: %d\n", neq1 - neq2);
    
    // 再次添加边界条件
    RgLog("--- Applying BC: Fix node 1 ---\n");
    mesh.Node(1).SetDofState(uIdx, DOF_PRESCRIBED);
    mesh.Node(1).SetDofState(vIdx, DOF_PRESCRIBED);
    mesh.Node(1).SetDofState(wIdx, DOF_PRESCRIBED);
    
    // 再次重新编号
    RgLog("--- Renumbering again ---\n");
    int neq3 = numbering.RenumberEquations();
    RgLog("Equations: %d\n\n", neq3);
    
    numbering.Validate();
}

//=============================================================================
// 示例5: 查询方程号
//=============================================================================

void Example5_QueryEquations()
{
    RgLog("\n");
    RgLog("=========================================================\n");
    RgLog("Example 5: Querying Equation Numbers\n");
    RgLog("=========================================================\n\n");
    
    FEModel model;
    FEMesh& mesh = model.GetMesh();
    
    mesh.CreateNodes(8);
    
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(1, ElementType::FE_HEX8G8);
    mesh.AddDomain(dom);
    
    model.Init();
    
    // 固定节点0
    RgDofSchema& schema = model.GetDofSchema();
    mesh.Node(0).FixDof(0, 0.0);
    mesh.Node(0).FixDof(1, 0.0);
    mesh.Node(0).FixDof(2, 0.0);
    
    // 编号
    FEEquationNumbering numbering(&model);
    numbering.NumberEquations();
    
    // 查询各节点的方程号
    RgLog("Node equation numbers:\n");
    RgLog("-------------------------------------------------------\n");
    
    int dofsPerNode = schema.GetDofsPerNode();
    
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        const FENode& node = mesh.Node(i);
        
        RgLog("Node %d: ", i);
        
        for (int d = 0; d < dofsPerNode; ++d)
        {
            int eqn = node.GetEquationNumber(d);
            
            DofType dofType = schema.GetDofType(d);
            std::string dofName = schema.GetDofName(dofType);
            
            if (eqn >= 0)
            {
                RgLog("%s=%d ", dofName.c_str(), eqn);
            }
            else
            {
                RgLog("%s=- ", dofName.c_str());
            }
        }
        
        RgLog("\n");
    }
    
    RgLog("-------------------------------------------------------\n\n");
    
    // 使用编号器查询
    RgLog("Using FEEquationNumbering::GetEquation():\n");
    
    int nodeId = 5;
    int dofIdx = schema.GetDofIndex("u");
    
    int eqn = numbering.GetEquation(nodeId, dofIdx);
    
    RgLog("Node %d, DOF 'u' (index %d) -> Equation %d\n\n",
          nodeId, dofIdx, eqn);
}

//=============================================================================
// 主函数：运行所有示例
//=============================================================================

int main()
{
    RgLog("======================================================\n");
    RgLog("FEEquationNumbering Usage Examples\n");
    RgLog("======================================================\n");
    
    // 运行示例
    Example1_BasicSolidModel();
    Example2_MixedSolidShellModel();
    Example3_ComplexBoundaryConditions();
    Example4_Renumbering();
    Example5_QueryEquations();
    
    RgLog("======================================================\n");
    RgLog("All examples completed\n");
    RgLog("======================================================\n");
    
    return 0;
}
