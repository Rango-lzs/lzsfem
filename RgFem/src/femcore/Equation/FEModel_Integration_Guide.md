/*********************************************************************
 * \file   FEModel_Integration_Guide.md
 * \brief  How to integrate equation numbering into FEModel
 *
 * \date   February 2026
 *********************************************************************/

# 方程编号集成到FEModel的完整指南

## 修改FEModel类

### 1. 在FEModel.h中添加

```cpp
#include "FEEquationNumbering.h"

class FEModel
{
public:
    // ... 现有成员 ...
    
    /// 获取方程编号器
    FEEquationNumbering& GetEquationNumbering() { return m_eqnNumbering; }
    
    /// 获取总方程数
    int GetTotalEquations() const { return m_eqnNumbering.GetTotalEquations(); }
    
    /// 执行方程编号
    int NumberEquations();
    
    /// 重新编号（边界条件改变后）
    int RenumberEquations();
    
private:
    RgDofSchema m_dofSchema;                // DOF配置
    FEEquationNumbering m_eqnNumbering;     // 方程编号器（新增）
    
    // ... 其他成员 ...
};
```

### 2. 在FEModel.cpp中实现

```cpp
FEModel::FEModel()
    : m_dofSchema()
    , m_eqnNumbering(this)  // 传入this指针
{
    // ... 其他初始化 ...
}

int FEModel::NumberEquations()
{
    return m_eqnNumbering.NumberEquations();
}

int FEModel::RenumberEquations()
{
    return m_eqnNumbering.RenumberEquations();
}
```

### 3. 修改FEModel::Init()

```cpp
bool FEModel::Init()
{
    RgLog("\n=== Initializing FE Model ===\n");
    
    FEMesh& mesh = GetMesh();
    
    //=========================================================================
    // STEP 1: 配置DOF Schema
    //=========================================================================
    
    RgLog("Step 1: Configuring DOF schema...\n");
    
    if (!AutoConfigureDofSchema())
    {
        RgLogError("Failed to configure DOF schema");
        return false;
    }
    
    //=========================================================================
    // STEP 2: 为所有节点分配DOF空间
    //=========================================================================
    
    RgLog("Step 2: Allocating node DOFs...\n");
    
    int dofsPerNode = m_dofSchema.GetDofsPerNode();
    
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        node.SetDOFS(dofsPerNode);
    }
    
    //=========================================================================
    // STEP 3: 激活节点DOF（根据单元连接性）
    //=========================================================================
    
    RgLog("Step 3: Activating node DOFs...\n");
    
    ActivateNodeDofs();
    
    //=========================================================================
    // STEP 4: 初始化材料、域等
    //=========================================================================
    
    RgLog("Step 4: Initializing materials...\n");
    
    if (!InitMaterials())
    {
        RgLogError("Failed to initialize materials");
        return false;
    }
    
    if (!InitMesh())
    {
        RgLogError("Failed to initialize mesh");
        return false;
    }
    
    //=========================================================================
    // STEP 5: 应用初始条件和边界条件
    //=========================================================================
    
    RgLog("Step 5: Applying initial and boundary conditions...\n");
    
    // 初始条件
    for (int i = 0; i < InitialConditions(); ++i)
    {
        if (!InitialCondition(i)->Init())
        {
            RgLogError("Failed to initialize IC %d", i);
            return false;
        }
    }
    
    // 边界条件
    if (!InitBCs())
    {
        RgLogError("Failed to initialize BCs");
        return false;
    }
    
    //=========================================================================
    // STEP 6: 全局方程编号 ⭐ 关键步骤
    //=========================================================================
    
    RgLog("Step 6: Numbering equations...\n");
    
    int neq = NumberEquations();
    
    if (neq <= 0)
    {
        RgLogError("No equations to solve!");
        return false;
    }
    
    RgLog("Total equations: %d\n", neq);
    
    //=========================================================================
    // 继续其他初始化
    //=========================================================================
    
    // 初始化接触
    if (!InitContact())
    {
        RgLogError("Failed to initialize contact");
        return false;
    }
    
    // 初始化载荷
    if (!InitModelLoads())
    {
        RgLogError("Failed to initialize model loads");
        return false;
    }
    
    RgLog("=== Model initialization complete ===\n\n");
    
    return true;
}
```

## 在求解器中使用

### 1. 创建全局矩阵和向量

```cpp
class FESolver
{
public:
    bool BuildMatrixProfile(FEGlobalMatrix& K)
    {
        FEModel* model = GetFEModel();
        
        // 获取总方程数
        int neq = model->GetTotalEquations();
        
        // 调整矩阵大小
        K.Create(neq, neq);
        
        // ... 构建矩阵profile ...
        
        return true;
    }
    
    bool Solve()
    {
        FEModel* model = GetFEModel();
        int neq = model->GetTotalEquations();
        
        // 创建求解向量
        std::vector<double> R(neq, 0.0);  // 残差向量
        std::vector<double> u(neq, 0.0);  // 解向量
        
        // 装配全局矩阵和向量
        // ...
        
        // 求解
        // ...
        
        // 更新节点DOF值
        UpdateSolution(u);
        
        return true;
    }
    
    void UpdateSolution(const std::vector<double>& u)
    {
        FEModel* model = GetFEModel();
        FEMesh& mesh = model->GetMesh();
        RgDofSchema& schema = model->GetDofSchema();
        
        int nNodes = mesh.Nodes();
        int dofsPerNode = schema.GetDofsPerNode();
        
        // 将解向量分配到节点
        for (int i = 0; i < nNodes; ++i)
        {
            FENode& node = mesh.Node(i);
            
            for (int d = 0; d < dofsPerNode; ++d)
            {
                int eqn = node.GetEquationNumber(d);
                
                if (eqn >= 0)  // 此DOF有方程号
                {
                    node.set(d, u[eqn]);
                }
                // 固定DOF保持原值
            }
        }
    }
};
```

### 2. 单元刚度装配

```cpp
void RgSolidElement::AssembleStiffness(FEGlobalMatrix& K)
{
    FEModel* model = GetFEModel();
    FEMesh& mesh = model->GetMesh();
    RgDofSchema& schema = model->GetDofSchema();
    
    int nNodes = NodeSize();
    
    // 获取solid单元需要的DOF索引
    RgDofSet solidDofs = RgDofSet::Solid3D();
    std::vector<int> dofIdx = schema.GetIndices(solidDofs);
    // dofIdx = {0, 1, 2} for u,v,w
    
    int elemDofsPerNode = dofIdx.size();  // 3
    
    // 计算单元刚度矩阵
    int totalElemDofs = nNodes * elemDofsPerNode;  // 24 for hex8
    Matrix Ke(totalElemDofs, totalElemDofs);
    ComputeElementStiffness(Ke);
    
    // 装配到全局矩阵
    for (int i = 0; i < nNodes; ++i)
    {
        int nodeId = getNodeId(i);
        FENode& node_i = mesh.Node(nodeId);
        
        for (int di = 0; di < elemDofsPerNode; ++di)
        {
            int dofIndex_i = dofIdx[di];
            int eqn_i = node_i.GetEquationNumber(dofIndex_i);
            
            // ⭐ 跳过没有方程号的DOF（固定或未激活）
            if (eqn_i < 0) continue;
            
            int localRow = i * elemDofsPerNode + di;
            
            for (int j = 0; j < nNodes; ++j)
            {
                int nodeId_j = getNodeId(j);
                FENode& node_j = mesh.Node(nodeId_j);
                
                for (int dj = 0; dj < elemDofsPerNode; ++dj)
                {
                    int dofIndex_j = dofIdx[dj];
                    int eqn_j = node_j.GetEquationNumber(dofIndex_j);
                    
                    // ⭐ 跳过没有方程号的DOF
                    if (eqn_j < 0) continue;
                    
                    int localCol = j * elemDofsPerNode + dj;
                    
                    // ⭐ 使用方程号装配
                    K.add(eqn_i, eqn_j, Ke(localRow, localCol));
                }
            }
        }
    }
}
```

### 3. 载荷向量装配

```cpp
void RgNodalLoad::AssembleLoad(std::vector<double>& R)
{
    FEModel* model = GetFEModel();
    FEMesh& mesh = model->GetMesh();
    RgDofSchema& schema = model->GetDofSchema();
    
    FENodeSet* nodeSet = GetNodeSet();
    
    int dofIdx = GetDofIndex();  // 施加载荷的DOF索引
    double loadValue = GetMagnitude();
    
    for (int i = 0; i < nodeSet->Size(); ++i)
    {
        int nodeId = nodeSet->Node(i);
        FENode& node = mesh.Node(nodeId);
        
        int eqn = node.GetEquationNumber(dofIdx);
        
        // ⭐ 只对有方程号的DOF施加载荷
        if (eqn >= 0)
        {
            R[eqn] += loadValue;
        }
    }
}
```

## 边界条件改变后重新编号

```cpp
bool FEModel::ApplyNewBoundaryCondition(int nodeId, int dofIndex, double value)
{
    FEMesh& mesh = GetMesh();
    FENode& node = mesh.Node(nodeId);
    
    // 应用边界条件
    node.SetDofState(dofIndex, DOF_PRESCRIBED);
    node.set(dofIndex, value);
    
    // ⭐ 重新编号
    RgLog("Boundary condition changed, renumbering...\n");
    int neq = RenumberEquations();
    
    RgLog("New equation count: %d\n", neq);
    
    return true;
}
```

## 诊断和调试

### 打印方程编号信息

```cpp
void FEModel::PrintEquationInfo()
{
    m_eqnNumbering.PrintStatistics();
    m_eqnNumbering.PrintDofDistribution();
}

void FEModel::ExportEquationInfo(const char* filename)
{
    m_eqnNumbering.ExportNumbering(filename);
}

void FEModel::ValidateEquations()
{
    if (!m_eqnNumbering.Validate())
    {
        RgLogError("Equation numbering validation failed!");
    }
}
```

## 完整的工作流程

```cpp
int main()
{
    // 1. 创建模型
    FEModel model;
    model.SetActiveModule("solid");
    
    FEMesh& mesh = model.GetMesh();
    
    // 2. 创建网格
    mesh.CreateNodes(100);
    // ... 设置节点坐标 ...
    
    // 3. 创建单元
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(10, ElementType::FE_HEX8G8);
    // ... 设置单元连接性 ...
    mesh.AddDomain(dom);
    
    // 4. 创建材料
    FEMaterial* mat = new FELinearElasticMaterial();
    model.AddMaterial(mat);
    dom->SetMaterial(mat);
    
    // 5. 应用边界条件（在Init之前或之后都可以）
    // 方式1: 直接操作节点
    RgDofSchema& schema = model.GetDofSchema();
    // ... 需要等Init之后才有schema ...
    
    // 6. 初始化模型 ⭐ 会自动执行方程编号
    if (!model.Init())
    {
        RgLogError("Model initialization failed");
        return -1;
    }
    
    // 7. 应用边界条件（在Init之后）
    FENodeSet* fixedNodes = mesh.FindNodeSet("FixedNodes");
    if (fixedNodes)
    {
        int uIdx = schema.GetDofIndex("u");
        int vIdx = schema.GetDofIndex("v");
        int wIdx = schema.GetDofIndex("w");
        
        for (int i = 0; i < fixedNodes->Size(); ++i)
        {
            int nodeId = fixedNodes->Node(i);
            FENode& node = mesh.Node(nodeId);
            
            node.SetDofState(uIdx, DOF_PRESCRIBED);
            node.SetDofState(vIdx, DOF_PRESCRIBED);
            node.SetDofState(wIdx, DOF_PRESCRIBED);
        }
        
        // ⭐ 重新编号
        model.RenumberEquations();
    }
    
    // 8. 打印方程信息
    model.PrintEquationInfo();
    model.ExportEquationInfo("equations.txt");
    
    // 9. 求解
    if (!model.Solve())
    {
        RgLogError("Solution failed");
        return -1;
    }
    
    return 0;
}
```

## 关键点总结

### ✅ 方程编号时机

1. **Init()结束时自动编号** - 推荐
2. **边界条件改变后重新编号** - 必须
3. **每个求解步开始前检查** - 可选

### ✅ 使用方程号的地方

1. **创建全局矩阵** - 使用neq确定大小
2. **装配刚度矩阵** - 使用eqn作为矩阵索引
3. **装配载荷向量** - 使用eqn作为向量索引
4. **求解后更新** - 使用eqn分配解到节点

### ✅ 重要检查

```cpp
int eqn = node.GetEquationNumber(dofIdx);
if (eqn < 0)
{
    // 此DOF没有方程号（固定或未激活）
    // 跳过，不要装配
    continue;
}
// 使用eqn装配
K.add(eqn, eqn, value);
```

### ❌ 常见错误

1. ❌ 忘记检查 `eqn < 0`
2. ❌ 边界条件改变后忘记重新编号
3. ❌ 使用 `nodeId * dofsPerNode + dofIdx` 代替方程号
4. ❌ 在Init()之前访问方程号

## 这就是完整的集成方案！
