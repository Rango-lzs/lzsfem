# RgAnalysisStep 系统设计文档

## 1. 概述

RgAnalysisStep 是一个全新的分析步骤管理系统，专门用于处理有限元分析中的多步骤问题。它解决了以下核心问题：

- **步骤间继承**: 边界条件和载荷可以从前一个步骤继承
- **步骤独立性**: 每个步骤可以定义自己的边界条件和载荷
- **灵活管理**: 支持添加、删除、修改步骤中的边界条件和载荷

## 2. 类层次结构

```
RgAnalysisStep (抽象基类)
    ├── RgStaticAnalysisStep (静态分析)
    └── RgDynamicAnalysisStep (动态分析)
```

## 3. 核心概念

### 3.1 边界条件和载荷的归属

每个 `RgAnalysisStep` 维护两个列表：

1. **本步骤定义的 BCs/Loads** (`m_boundaryConditions`, `m_loads`)
   - 在此步骤中新定义的边界条件和载荷
   - 由此步骤"拥有"（从管理角度）

2. **继承的 BCs/Loads** (`m_inheritedBCs`, `m_inheritedLoads`)
   - 从前面步骤继承的边界条件和载荷
   - 只是引用，不拥有所有权

### 3.2 激活模式

```cpp
enum class StepActivationMode
{
    NEW,         // 仅使用本步骤定义的 BCs/Loads
    INHERITED,   // 继承 + 本步骤定义 (默认)
    REPLACE      // 替换所有（只用本步骤定义）
};
```

**INHERITED 模式**（最常用）:
```
Step 1: BC1, Load1
Step 2: BC2, Load2
    → 实际激活: BC1, Load1, BC2, Load2

Step 3: BC3
    → 实际激活: BC1, Load1, BC2, Load2, BC3
```

**NEW 模式**:
```
Step 1: BC1, Load1
Step 2 (NEW): BC2, Load2
    → 实际激活: BC2, Load2 (不继承)
```

**REPLACE 模式**:
```
Step 1: BC1, Load1
Step 2 (REPLACE): BC2, Load2
    → 实际激活: BC2, Load2 (替换全部)
```

### 3.3 继承机制

`InheritFromPreviousStep()` 方法的工作流程:

1. 获取前一步骤的所有激活 BCs/Loads
2. 检查当前步骤是否重新定义了同名 BC/Load
3. 如果没有重新定义，则添加到继承列表
4. 如果重新定义了，则使用新定义（覆盖）

**示例**:
```cpp
// Step 1
step1->AddBoundaryCondition(bc_fixed);  // 名为 "FixedSupport"
step1->AddLoad(load_pressure);          // 名为 "Pressure"

// Step 2
step2->SetPreviousStep(step1);
step2->AddBoundaryCondition(bc_new);    // 名为 "NewConstraint"
step2->InheritFromPreviousStep();

// 结果: Step 2 激活的 BCs = [FixedSupport(继承), NewConstraint(新定义)]
//      Step 2 激活的 Loads = [Pressure(继承)]
```

## 4. 主要接口

### 4.1 RgAnalysisStep 基类

#### 基本属性
```cpp
void SetName(const std::string& name);
std::string GetName() const;
void SetStepNumber(int n);
int GetStepNumber() const;
```

#### 时间控制
```cpp
void SetTimePeriod(double t);              // 总时间周期
void SetInitialTimeIncrement(double dt);   // 初始时间增量
void SetMinTimeIncrement(double dt);       // 最小时间增量
void SetMaxTimeIncrement(double dt);       // 最大时间增量
```

#### 边界条件管理
```cpp
void AddBoundaryCondition(FEBoundaryCondition* bc);
void RemoveBoundaryCondition(const std::string& name);
FEBoundaryCondition* FindBoundaryCondition(const std::string& name);
int BoundaryConditions() const;  // 仅本步骤定义的数量
```

#### 载荷管理
```cpp
void AddLoad(FEModelLoad* load);
void RemoveLoad(const std::string& name);
FEModelLoad* FindLoad(const std::string& name);
int Loads() const;  // 仅本步骤定义的数量
```

#### 继承管理
```cpp
void SetActivationMode(StepActivationMode mode);
void SetPreviousStep(RgAnalysisStep* prevStep);
void InheritFromPreviousStep();

// 获取所有激活的 BCs/Loads（包括继承的）
std::vector<FEBoundaryCondition*> GetAllActiveBCs();
std::vector<FEModelLoad*> GetAllActiveLoads();
```

#### 生命周期
```cpp
virtual bool Initialize();   // 初始化步骤（包括继承处理）
virtual bool Activate();     // 激活所有 BCs 和 Loads
virtual bool Deactivate();   // 去激活本步骤定义的 BCs 和 Loads
virtual bool Solve() = 0;    // 求解（纯虚函数）
```

### 4.2 RgStaticAnalysisStep

```cpp
void SetNonlinear(bool b);
void SetMaxIterations(int n);
void SetConvergenceTolerance(double tol);
```

### 4.3 RgDynamicAnalysisStep

```cpp
void SetExplicit(bool b);
void SetAlpha(double alpha);    // Newmark 参数
void SetBeta(double beta);
void SetGamma(double gamma);
```

## 5. Abaqus INP 解析集成

### 5.1 parseStep 工作流程

```
parseStep()
    ↓
解析 *STEP 关键字参数
    ↓
创建 RgAnalysisStep 对象（根据 *STATIC/*DYNAMIC）
    ↓
循环解析步骤内容:
    *STATIC/*DYNAMIC → 设置时间参数
    *BOUNDARY → 调用 parseBoundary()，收集新 BCs
    *CLOAD → 调用 parseCload()，收集新 Loads
    *DLOAD → 调用 parseDload()，收集新 Loads
    *END STEP → 退出
    ↓
将收集的 BCs/Loads 添加到 step
    ↓
设置继承关系（如果有前一步骤）
    ↓
添加到 FEModel
```

### 5.2 BC/Load 收集机制

```cpp
// 存储解析前的 BC 数量
int bcCountBefore = m_model->BoundaryConditions();

// 解析 BOUNDARY
parseBoundary(file);  // 这会将新 BC 添加到 m_model

// 收集新添加的 BCs
int bcCountAfter = m_model->BoundaryConditions();
for (int i = bcCountBefore; i < bcCountAfter; ++i)
{
    FEBoundaryCondition* bc = m_model->BoundaryCondition(i);
    analysisStep->AddBoundaryCondition(bc);
}
```

### 5.3 Abaqus INP 示例

#### 示例 1: 两步加载

```inp
*Heading
Two-step loading analysis

*Node
1, 0.0, 0.0, 0.0
...

*Element, type=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
...

*Nset, nset=FixedNodes
1, 2, 3

*Nset, nset=LoadNodes
100, 101, 102

** ====================================
** STEP 1: Apply half load
** ====================================
*Step, name=HalfLoad
*Static
0.1, 1.0, 1e-5, 0.2

*Boundary
FixedNodes, 1, 3, 0.0

*Cload
LoadNodes, 2, -500.0

*End Step

** ====================================
** STEP 2: Apply full load
** ====================================
*Step, name=FullLoad
*Static
0.1, 1.0, 1e-5, 0.2

** Boundary from Step 1 is inherited
** Add additional load
*Cload
LoadNodes, 2, -1000.0

*End Step
```

**解析结果**:

```
Step 1 "HalfLoad":
  - New BCs: FixedNodes (U1=U2=U3=0)
  - New Loads: LoadNodes (F2=-500)
  - Active BCs: [FixedNodes]
  - Active Loads: [LoadNodes F2=-500]

Step 2 "FullLoad":
  - New BCs: (none)
  - New Loads: LoadNodes (F2=-1000, 覆盖前一个)
  - Inherited BCs: [FixedNodes]
  - Inherited Loads: (none, 被覆盖)
  - Active BCs: [FixedNodes]
  - Active Loads: [LoadNodes F2=-1000]
```

#### 示例 2: 多个边界条件步骤

```inp
*Step, name=Prestress
*Static
0.05, 0.5

*Boundary
FixedSupport, 1, 3, 0.0

*Cload
PreloadNodes, 3, -100.0

*End Step

*Step, name=MainLoad
*Static
0.1, 1.0

** Inherit FixedSupport and PreloadNodes
** Add new boundary condition
*Boundary
RollerSupport, 1, 1, 0.0
RollerSupport, 3, 3, 0.0

** Add main load
*Cload
MainLoadNodes, 2, -1000.0

*End Step

*Step, name=Unload
*Static
0.1, 0.5

** Remove PreloadNodes by setting to 0
*Cload
PreloadNodes, 3, 0.0

** Keep other BCs and loads
*End Step
```

**解析结果**:

```
Step 1 "Prestress":
  - Active BCs: [FixedSupport]
  - Active Loads: [PreloadNodes F3=-100]

Step 2 "MainLoad":
  - Inherited: FixedSupport, PreloadNodes
  - New: RollerSupport, MainLoadNodes
  - Active BCs: [FixedSupport, RollerSupport]
  - Active Loads: [PreloadNodes F3=-100, MainLoadNodes F2=-1000]

Step 3 "Unload":
  - Inherited: FixedSupport, RollerSupport, MainLoadNodes
  - Modified: PreloadNodes F3=0 (实际是移除)
  - Active BCs: [FixedSupport, RollerSupport]
  - Active Loads: [MainLoadNodes F2=-1000]
```

## 6. 使用示例

### 6.1 手动创建步骤

```cpp
FEModel model;

// 创建第一个步骤
RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
step1->SetName("InitialLoad");
step1->SetTimePeriod(1.0);
step1->SetInitialTimeIncrement(0.1);

// 添加边界条件
FEBoundaryCondition* bc_fixed = ...; // 创建固定支撑
step1->AddBoundaryCondition(bc_fixed);

// 添加载荷
FEModelLoad* load1 = ...; // 创建载荷
step1->AddLoad(load1);

// 添加到模型
model.AddStep(step1);

// 创建第二个步骤
RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
step2->SetName("AdditionalLoad");
step2->SetTimePeriod(1.0);
step2->SetPreviousStep(step1);  // 设置前一步骤
step2->SetActivationMode(StepActivationMode::INHERITED);

// 添加新载荷
FEModelLoad* load2 = ...; // 创建额外载荷
step2->AddLoad(load2);

// 添加到模型
model.AddStep(step2);

// 初始化和求解
model.Init();  // 会调用每个 step 的 Initialize()
model.Solve(); // 会依次调用每个 step 的 Solve()
```

### 6.2 从 Abaqus 文件加载

```cpp
FEModel model;
AbaqusImport importer;

// 加载 INP 文件（会自动创建步骤）
if (!importer.load("model.inp", &model))
{
    RgLogError("Failed to load INP file\n");
    return false;
}

// 检查步骤
int nsteps = model.Steps();
RgLog("Loaded %d steps\n", nsteps);

for (int i = 0; i < nsteps; ++i)
{
    RgAnalysisStep* step = dynamic_cast<RgAnalysisStep*>(model.GetStep(i));
    if (step)
    {
        RgLog("Step %d: %s\n", i+1, step->GetName().c_str());
        
        // 获取所有激活的 BCs
        auto activeBCs = step->GetAllActiveBCs();
        RgLog("  Active BCs: %d\n", (int)activeBCs.size());
        
        // 获取所有激活的载荷
        auto activeLoads = step->GetAllActiveLoads();
        RgLog("  Active loads: %d\n", (int)activeLoads.size());
    }
}

// 初始化并求解
model.Init();
model.Solve();
```

### 6.3 步骤间修改

```cpp
// 获取某个步骤
RgAnalysisStep* step = dynamic_cast<RgAnalysisStep*>(model.GetStep(1));

// 添加新的边界条件
FEBoundaryCondition* newBC = ...;
step->AddBoundaryCondition(newBC);

// 移除某个载荷
step->RemoveLoad("OldLoadName");

// 查找并修改边界条件
FEBoundaryCondition* bc = step->FindBoundaryCondition("FixedSupport");
if (bc)
{
    // 修改 BC 参数...
}

// 重新初始化继承
step->InheritFromPreviousStep();
```

## 7. 高级特性

### 7.1 条件继承

可以重写 `InheritFromPreviousStep()` 来实现自定义继承逻辑:

```cpp
class MyCustomStep : public RgStaticAnalysisStep
{
public:
    void InheritFromPreviousStep() override
    {
        // 只继承名称包含 "Permanent" 的 BCs
        if (m_previousStep)
        {
            auto prevBCs = m_previousStep->GetAllActiveBCs();
            for (auto bc : prevBCs)
            {
                if (bc->GetName().find("Permanent") != std::string::npos)
                {
                    if (!FindBoundaryCondition(bc->GetName()))
                    {
                        m_inheritedBCs.push_back(bc);
                    }
                }
            }
        }
    }
};
```

### 7.2 步骤组

可以创建步骤组来管理相关步骤:

```cpp
class RgStepGroup
{
    std::vector<RgAnalysisStep*> m_steps;
    
public:
    void AddStep(RgAnalysisStep* step) 
    { 
        m_steps.push_back(step); 
    }
    
    bool InitializeAll()
    {
        for (auto step : m_steps)
        {
            if (!step->Initialize())
                return false;
        }
        return true;
    }
    
    bool SolveAll()
    {
        for (auto step : m_steps)
        {
            if (!step->Solve())
                return false;
        }
        return true;
    }
};
```

## 8. 调试和日志

### 8.1 详细日志输出

当启用详细日志时，`Initialize()` 会输出:

```
Initializing step 'MainLoad'
  Inherited BCs: 2
  New BCs: 1
  Total active BCs: 3
  Inherited loads: 1
  New loads: 2
  Total active loads: 3
```

### 8.2 验证步骤配置

```cpp
void ValidateStep(RgAnalysisStep* step)
{
    RgLog("Validating step '%s'\n", step->GetName().c_str());
    
    // 检查时间参数
    if (step->GetTimePeriod() <= 0)
    {
        RgLogWarning("  Time period <= 0\n");
    }
    
    if (step->GetInitialTimeIncrement() <= 0)
    {
        RgLogWarning("  Initial time increment <= 0\n");
    }
    
    // 检查是否有 BCs 或载荷
    auto bcs = step->GetAllActiveBCs();
    auto loads = step->GetAllActiveLoads();
    
    if (bcs.empty())
    {
        RgLogWarning("  No boundary conditions!\n");
    }
    
    if (loads.empty())
    {
        RgLogWarning("  No loads!\n");
    }
    
    RgLog("  Validation complete\n");
}
```

## 9. 性能优化建议

1. **避免重复继承**: 只在 `Initialize()` 中调用一次 `InheritFromPreviousStep()`

2. **使用引用而非拷贝**: 继承列表存储指针，不拷贝对象

3. **延迟激活**: 只在 `Activate()` 时激活 BCs/Loads，而不是在添加时

4. **批量操作**: 如果要添加多个 BCs，先全部添加，最后统一继承

## 10. 注意事项

### 10.1 内存管理

- `RgAnalysisStep` **不拥有** BC 和 Load 对象
- BC 和 Load 由 `FEModel` 管理和释放
- 继承列表只是引用，析构时不要删除

### 10.2 步骤顺序

- 步骤必须按顺序添加到 FEModel
- 前一步骤必须先于当前步骤添加
- `SetPreviousStep()` 应在 `InheritFromPreviousStep()` 之前调用

### 10.3 名称冲突

- 如果当前步骤定义了与继承 BC 同名的 BC，使用当前步骤的定义
- 这实现了"覆盖"语义

## 11. 未来扩展

可能的扩展方向:

1. **步骤依赖**: 步骤可以依赖多个前序步骤
2. **条件激活**: 根据前一步骤的结果决定是否激活
3. **子步骤**: 支持步骤嵌套
4. **步骤模板**: 预定义常用步骤配置
5. **并行步骤**: 支持独立步骤的并行求解

## 12. 总结

RgAnalysisStep 系统提供了:

✅ **清晰的归属关系**: 每个 BC/Load 明确属于某个步骤  
✅ **灵活的继承机制**: 支持多种继承模式  
✅ **与 Abaqus 兼容**: 完整支持 Abaqus INP 步骤定义  
✅ **易于扩展**: 基于继承的设计便于添加新类型  
✅ **健壮的错误处理**: 详细的日志和验证机制  

这个设计为多步骤有限元分析提供了坚实的基础。
