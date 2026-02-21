# StaticAnalysis 开发计划

## 项目现状分析

### 已完成的框架

| 模块 | 状态 | 说明 |
|------|------|------|
| **FEModel / FEMesh / FENode** | ✅ 基本完成 | 模型、网格、节点管理已实现 |
| **AnalysisStep / StaticStep** | ✅ 流程完成 | 线性/非线性分析步骤控制逻辑已实现 |
| **RgAnalysis** | ✅ 框架完成 | 分析管理器，支持创建 StaticAnalysis |
| **FENewtonSolver** | ✅ 迭代逻辑完成 | Newton-Raphson 迭代、收敛检查、线搜索 |
| **FEStaticSolver** | ⚠️ 骨架完成 | 关键装配调用被注释掉 |
| **LinearSolver** | ✅ 多个可用 | SkylineSolver、PardisoSolver、LUSolver |
| **RgSolidDomain** | ⚠️ 接口完成 | `StiffnessMatrix()`/`InternalForces()` 委托给 Assembler |
| **RgElasticAssembler** | ⚠️ 框架完成 | 装配循环已写，但依赖元素级计算 |
| **Element (2D/3D)** | ❌ 多为占位 | 大部分 `calculateStiffnessMatrix()` 是 placeholder |
| **Material** | ⚠️ 部分实现 | 线弹性材料模型需确认是否完整 |
| **BoundaryCondition / Load** | ⚠️ 接口完成 | 应用到方程系统的连接需验证 |
| **Equation Numbering** | ✅ 已实现 | `RgDofSchema` + `FEEquationNumbering` |
| **Input (Abaqus/FEBio)** | ⚠️ 部分实现 | 可读取基本网格，材料/BC转换不完整 |

### 关键断点（数据流中断的地方）

```
RgAnalysis::run()
  → StaticStep::execute()
    → FEStaticSolver::SolveStep()  [调用 Quasin()]
      → FEStaticSolver::StiffnessMatrix()
        → dom.StiffnessMatrix(LS)     ← ❌ 被注释掉！
          → Assembler::ElementStiffness()
            → Element::calculateStiffnessMatrix()  ← ❌ 多为空实现！
      → FEStaticSolver::Residual()
        → InternalForces()            ← ❌ 被注释掉！
        → ExternalForces()            ← ❌ 被注释掉！
      → LinearSolver::Factor() + BackSolve()  ← ✅ 可用
```

---

## 开发计划（分阶段）

### 阶段 0：选择最简单的验证问题
**目标**：定义一个可以手算验证的 benchmark  
**预计工作量**：0.5 天

选用 **3D 悬臂梁单元素测试**（单个 Hex8 单元）：
- 1 个 Hex8 单元，8 个节点
- 线弹性材料 (E=200GPa, ν=0.3)
- 一端固定（4个节点全约束）
- 另一端施加集中力
- 解析解可手算验证

```
固定端           自由端
4-----7          施加力 F
|     |    →     方向：Z
5-----6
|     |
0-----3
|     |
1-----2
```

---

### 阶段 1：打通 Element 级计算
**目标**：让一个 3D 实体单元能正确计算刚度矩阵和内力  
**预计工作量**：3-5 天

#### 1.1 实现 Hex8 单元刚度矩阵
**文件**：选择在 `RgElasticAssembler` 或对应的 3D 实体元素类中实现

```cpp
// 需要实现的核心函数
void ElementStiffness(int iel, FEElementMatrix& ke)
{
    // 1. 获取单元节点坐标
    // 2. 8 点 Gauss 积分
    // 3. 对每个积分点：
    //    a. 计算形函数及其导数 (N, dN/dξ)
    //    b. 计算 Jacobian J = dN/dξ * x_node
    //    c. 计算 B 矩阵 (6×24 应变-位移矩阵)
    //    d. 获取材料本构矩阵 D (6×6)
    //    e. Ke += B^T * D * B * det(J) * w_gauss
}
```

**检查清单**：
- [ ] `RgSolidElement` 能正确返回节点坐标
- [ ] Hex8 的形函数和导数计算正确（`FEElementTraits` 中可能已有）
- [ ] Jacobian 计算与逆变换
- [ ] B 矩阵组装（线性 = 小变形假设）
- [ ] 材料 D 矩阵接口（从 `RgMaterial` 获取）

#### 1.2 实现线弹性本构矩阵
**文件**：材料类中

```cpp
// 各向同性线弹性 D 矩阵 (6×6)
void GetElasticityTensor(double E, double nu, Matrix& D)
{
    double lambda = E * nu / ((1+nu) * (1-2*nu));
    double mu = E / (2*(1+nu));
    // 填充 D 矩阵
}
```

**检查清单**：
- [ ] 确认 `RgMaterial` 或 `FELinearElasticMaterial` 是否有工作的实现
- [ ] D 矩阵的 Voigt notation 约定与 B 矩阵一致

#### 1.3 实现单元内力向量
```cpp
void ElementInternalForces(RgSolidElement& el, vector<double>& fe)
{
    // fe = ∫ B^T * σ dV
    // 对线性问题：σ = D * ε = D * B * u_e
}
```

#### 1.4 单元级单元测试
编写独立测试，验证单个 Hex8 单元的：
- 刚度矩阵对称性
- 刚度矩阵秩 = 24 - 6 = 18（扣除刚体模态）
- patch test 通过

---

### 阶段 2：打通 Domain → Solver 装配链路
**目标**：让 FEStaticSolver 能正确装配全局刚度矩阵和载荷向量  
**预计工作量**：2-3 天

#### 2.1 取消注释并连接 Domain 调用
**文件**：`StaticSolver.cpp`

```cpp
// StiffnessMatrix() 中取消注释
for (int i = 0; i < mesh.Domains(); ++i)
{
    RgDomain& dom = mesh.Domain(i);
    if (dom.isActive())
    {
        dom.StiffnessMatrix(LS);  // ← 取消注释
    }
}
```

```cpp
// Residual() 中确保 InternalForces 和 ExternalForces 被正确调用
InternalForces(RHS);   // ← 确保内部有实际调用 domain 级别
ExternalForces(RHS);   // ← 确保外力正确装配
```

#### 2.2 确认 UnpackLM（方程号映射）正确
**文件**：`RgSolidDomain::UnpackLM()`

验证 `lm` 向量正确映射了：节点的 DOF → 全局方程号。
当前实现直接用 `node.getDofs()[j]`，需确认这返回的是全局方程号。

#### 2.3 确认 FEGlobalMatrix 装配正确
- `FELinearSystem::Assemble(ke)` 能正确将单元矩阵加到全局矩阵
- 全局矩阵的 profile（稀疏结构）正确建立

#### 2.4 确认 FEGlobalVector 装配正确
- 外力向量通过 `FEGlobalVector::Assemble()` 正确装配
- 节点载荷正确映射到方程号

---

### 阶段 3：打通 Solver 初始化和求解
**目标**：让线性求解器能正确初始化、分解、回代  
**预计工作量**：2-3 天

#### 3.1 LinearSolver 创建与关联
**文件**：`FENewtonSolver::Init()` 或 `FEStaticSolver::Init()`

```cpp
// 需要确保 m_plinsolve 被正确创建
// 当前代码中有 TODO 标记
if (m_plinsolve == nullptr)
{
    m_plinsolve = new SkylineSolver();  // 或 PardisoSolver
}
```

**检查清单**：
- [ ] `AllocateLinearSystem()` 正确创建 `FEGlobalMatrix` 和 `SparseMatrix`
- [ ] `CreateStiffness(true)` 正确建立矩阵 profile
- [ ] `ReformStiffness()` → `Factor()` → `BackSolve()` 链路畅通
- [ ] SkylineSolver 作为首选（不依赖外部库）

#### 3.2 方程编号系统
**文件**：`FESolver::InitEquations()`

确认：
- [ ] 自由 DOF 被分配正方程号（≥ 0）
- [ ] 固定 DOF 被分配负方程号（< 0）
- [ ] `m_neq` 正确反映总方程数

#### 3.3 边界条件应用
确认：
- [ ] 固定 BC：对应 DOF 被标记为 prescribed，方程号为负
- [ ] 集中力：正确加到载荷向量对应方程号位置
- [ ] prescribed displacement 的修正项（`m_Fd` 向量）正确处理

---

### 阶段 4：编写端到端测试
**目标**：一个完整的从模型建立到结果输出的测试程序  
**预计工作量**：1-2 天

#### 4.1 编写 TestStaticAnalysis 程序

```cpp
int main()
{
    // === 1. 创建模型 ===
    FEModel model;
    model.SetActiveModule("solid");
    FEMesh& mesh = model.GetMesh();
    
    // === 2. 创建节点（8个节点，1个Hex8单元）===
    mesh.CreateNodes(8);
    // 设置坐标：1×1×1 的立方体
    // Node 0: (0,0,0), Node 1: (1,0,0), ...
    
    // === 3. 创建单元 ===
    RgSolidDomain* dom = new RgSolidDomain(&model);
    dom->Create(1, ElementType::FE_HEX8G8);
    // 设置连接性 [0,1,2,3,4,5,6,7]
    mesh.AddDomain(dom);
    
    // === 4. 创建材料 ===
    // E = 200e3, nu = 0.3
    RgMaterial* mat = createLinearElasticMaterial(200e3, 0.3);
    model.AddMaterial(mat);
    dom->SetMaterial(mat);
    
    // === 5. 初始化 ===
    model.Init();
    
    // === 6. 应用边界条件 ===
    // 固定 z=0 面（节点 0,1,2,3）的所有DOF
    // 在 z=1 面（节点 4,5,6,7）施加 Z 方向集中力
    
    // === 7. 创建分析并运行 ===
    RgAnalysis analysis(&model);
    analysis.createStaticAnalysis("Test", 1.0);
    bool success = analysis.run();
    
    // === 8. 输出结果 ===
    // 打印各节点位移
    // 与解析解对比
    
    return success ? 0 : 1;
}
```

#### 4.2 验证标准
- 位移结果与手算/解析解误差 < 1%
- 反力平衡：∑F_reaction = ∑F_applied
- 应力均匀（单元素情况）

---

### 阶段 5：扩展验证
**目标**：多单元 patch test 和收敛性验证  
**预计工作量**：2-3 天

- [ ] 多单元悬臂梁（2×1×1, 4×1×1, 8×1×1）
- [ ] 验证 h-收敛（位移收敛到解析解）
- [ ] 应力 patch test（均匀应力场精确恢复）
- [ ] 从 Abaqus INP 文件读取简单模型并求解

---

## 优先级排序与依赖关系

```
阶段 0 (定义问题)
  │
  ▼
阶段 1 (单元级计算)  ← 核心工作量最大
  │
  ├─── 1.1 Hex8 刚度矩阵
  ├─── 1.2 线弹性 D 矩阵  
  ├─── 1.3 单元内力
  └─── 1.4 单元测试
  │
  ▼
阶段 2 (装配链路)
  │
  ├─── 2.1 取消注释 Domain 调用
  ├─── 2.2 验证 UnpackLM
  ├─── 2.3 全局矩阵装配
  └─── 2.4 全局载荷装配
  │
  ▼
阶段 3 (求解器初始化)
  │
  ├─── 3.1 LinearSolver 创建
  ├─── 3.2 方程编号
  └─── 3.3 边界条件
  │
  ▼
阶段 4 (端到端测试)
  │
  ▼
阶段 5 (扩展验证)
```

## 预计总工作量

| 阶段 | 工作量 | 累计 |
|------|--------|------|
| 阶段 0 | 0.5 天 | 0.5 天 |
| 阶段 1 | 3-5 天 | 4-6 天 |
| 阶段 2 | 2-3 天 | 6-9 天 |
| 阶段 3 | 2-3 天 | 8-12 天 |
| 阶段 4 | 1-2 天 | 9-14 天 |
| 阶段 5 | 2-3 天 | 11-17 天 |

**总计约 2-3 周**，核心突破在阶段 1（单元计算）和阶段 2（装配连接）。

## 建议的开发策略

1. **自底向上**：先确保单元级计算正确（可独立单元测试），再逐层向上连接
2. **最小路径优先**：只实现 Hex8 + 线弹性 + 集中力 + 固定BC，不要同时处理多种单元/材料/载荷
3. **频繁打印调试**：在关键节点打印矩阵大小、非零元数、范数等，快速定位问题
4. **利用已有框架**：项目中 `FEElementTraits` 已有 Hex8 的 Gauss 积分点和形函数数据，优先复用
