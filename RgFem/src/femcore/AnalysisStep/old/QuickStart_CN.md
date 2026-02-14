# RgAnalysisStep 快速入门指南

## 1. 文件清单

实现 RgAnalysisStep 系统需要以下文件:

```
RgAnalysisStep.h              - 类定义
RgAnalysisStep.cpp            - 类实现
AbaqusImport_parseStep_v2.cpp - 更新的 parseStep 实现
test_RgAnalysisStep.cpp       - 测试程序 (可选)
CMakeLists_RgAnalysisStep.txt - CMake 配置 (可选)
```

## 2. 集成步骤

### 步骤 1: 添加文件到项目

将以下文件添加到你的项目中:
- `RgAnalysisStep.h` → 放到 `femcore/FEAnalysis/` 或类似目录
- `RgAnalysisStep.cpp` → 同上

### 步骤 2: 更新 AbaqusImport

在 `AbaqusImport.cpp` 中:

1. 包含新头文件:
```cpp
#include "RgAnalysisStep.h"
```

2. 替换 `parseStep()` 函数:
   - 用 `AbaqusImport_parseStep_v2.cpp` 中的实现替换现有实现

3. 在 `AbaqusImport.h` 中添加:
```cpp
private:
    bool skipToNextKeyword(std::ifstream& file);
```

### 步骤 3: 更新 FEModel (如果需要)

如果你的 `FEModel::AddStep()` 接受 `FEAnalysis*` 参数，需要确保:

```cpp
// 在 FEModel.h 中
#include "RgAnalysisStep.h"

// FEModel::AddStep() 应该能接受 RgAnalysisStep*
void AddStep(FEAnalysis* step);  // 确保这个声明存在
```

### 步骤 4: 编译测试

```bash
# 添加到你的 CMakeLists.txt
add_library(RgAnalysisStep RgAnalysisStep.cpp)
target_link_libraries(YourProject RgAnalysisStep)

# 编译
cmake --build .
```

## 3. 基本使用

### 示例 1: 手动创建步骤

```cpp
#include "femcore/FEModel.h"
#include "RgAnalysisStep.h"

int main()
{
    // 创建模型
    FEModel model;
    model.SetActiveModule("solid");
    
    // 创建网格 (省略细节)
    // ...
    
    // 创建第一个分析步
    RgStaticAnalysisStep* step1 = new RgStaticAnalysisStep(&model);
    step1->SetName("ApplyLoad");
    step1->SetTimePeriod(1.0);
    step1->SetInitialTimeIncrement(0.1);
    
    // 添加边界条件 (假设已创建)
    step1->AddBoundaryCondition(bc_fixed);
    step1->AddLoad(load_pressure);
    
    // 添加到模型
    model.AddStep(step1);
    
    // 创建第二个步骤 (继承第一个步骤)
    RgStaticAnalysisStep* step2 = new RgStaticAnalysisStep(&model);
    step2->SetName("IncreaseLoad");
    step2->SetTimePeriod(1.0);
    step2->SetPreviousStep(step1);
    step2->SetActivationMode(StepActivationMode::INHERITED);
    
    // 添加额外载荷
    step2->AddLoad(load_additional);
    
    // 添加到模型
    model.AddStep(step2);
    
    // 初始化并求解
    model.Init();
    model.Solve();
    
    return 0;
}
```

### 示例 2: 从 Abaqus INP 加载

```cpp
#include "femcore/FEModel.h"
#include "AbaqusImport.h"
#include "RgAnalysisStep.h"

int main()
{
    // 创建模型
    FEModel model;
    
    // 加载 INP 文件
    AbaqusImport importer;
    if (!importer.load("model.inp", &model))
    {
        std::cerr << "Failed to load INP file\n";
        return 1;
    }
    
    // 检查加载的步骤
    int nsteps = model.Steps();
    std::cout << "Loaded " << nsteps << " steps\n";
    
    for (int i = 0; i < nsteps; ++i)
    {
        RgAnalysisStep* step = dynamic_cast<RgAnalysisStep*>(model.GetStep(i));
        if (step)
        {
            std::cout << "Step " << i+1 << ": " << step->GetName() << "\n";
            std::cout << "  Active BCs: " << step->GetAllActiveBCs().size() << "\n";
            std::cout << "  Active Loads: " << step->GetAllActiveLoads().size() << "\n";
        }
    }
    
    // 初始化并求解
    model.Init();
    model.Solve();
    
    return 0;
}
```

## 4. 测试验证

### 运行测试程序

```bash
# 编译测试程序
g++ -o test_step test_RgAnalysisStep.cpp RgAnalysisStep.cpp \
    -I/path/to/include -L/path/to/lib -lFEMCore

# 运行测试
./test_step
```

### 预期输出示例

```
=================================================
RgAnalysisStep Test Suite
=================================================

=================================================
Test 1: Basic Step Creation
=================================================
Created step: Step1-InitialLoad
  Time period: 1
  Initial dt: 0.1

Total steps in model: 2
Test 1 PASSED

...
```

## 5. 调试技巧

### 启用详细日志

确保你的日志系统启用:

```cpp
// 在初始化前
RgLog::SetVerbosity(VERBOSE);  // 根据你的日志系统调整
```

### 检查步骤配置

```cpp
void PrintStepInfo(RgAnalysisStep* step)
{
    std::cout << "Step: " << step->GetName() << "\n";
    std::cout << "  Number: " << step->GetStepNumber() << "\n";
    std::cout << "  Time period: " << step->GetTimePeriod() << "\n";
    std::cout << "  Initial dt: " << step->GetInitialTimeIncrement() << "\n";
    
    auto bcs = step->GetAllActiveBCs();
    std::cout << "  Active BCs: " << bcs.size() << "\n";
    for (auto bc : bcs)
    {
        std::cout << "    - " << bc->GetName() << "\n";
    }
    
    auto loads = step->GetAllActiveLoads();
    std::cout << "  Active Loads: " << loads.size() << "\n";
    for (auto load : loads)
    {
        std::cout << "    - " << load->GetName() << "\n";
    }
}
```

## 6. 常见问题

### Q1: 步骤没有继承前一个步骤的 BC/Load?

**A**: 检查:
1. 是否调用了 `SetPreviousStep()`
2. 是否设置了正确的 `ActivationMode`
3. 是否在 `Initialize()` 中调用了 `InheritFromPreviousStep()`

```cpp
step2->SetPreviousStep(step1);
step2->SetActivationMode(StepActivationMode::INHERITED);
step2->Initialize();  // 这会调用 InheritFromPreviousStep()
```

### Q2: Abaqus INP 解析后步骤中没有 BC/Load?

**A**: 检查:
1. `parseBoundary()`, `parseCload()`, `parseDload()` 是否正确实现
2. 这些函数是否将 BC/Load 添加到了 `FEModel`
3. `parseStep()` 中的收集逻辑是否正确

### Q3: 编译错误 - 找不到 RgAnalysisStep

**A**: 
1. 检查 `#include` 路径是否正确
2. 确保 `RgAnalysisStep.cpp` 已添加到编译列表
3. 检查 CMakeLists.txt 配置

### Q4: 运行时崩溃

**A**: 
1. 检查是否有空指针（`m_model`, `m_previousStep`）
2. 确保 BC/Load 对象在步骤生命周期内有效
3. 使用调试器检查调用栈

## 7. 下一步

1. **集成求解器**: 在 `RgStaticAnalysisStep::Solve()` 中实现实际求解逻辑

2. **扩展 BC/Load 类型**: 根据需要添加更多边界条件和载荷类型

3. **实现输出**: 在步骤中添加输出请求处理

4. **性能优化**: 对于大规模模型，优化继承机制

5. **添加验证**: 实现步骤配置的完整性检查

## 8. 获取帮助

如果遇到问题:

1. 查看详细文档: `RgAnalysisStep_Documentation_CN.md`
2. 运行测试程序: `test_RgAnalysisStep`
3. 检查日志输出
4. 使用调试器逐步执行

## 9. 许可和贡献

(根据你的项目添加许可信息)

---

**最后更新**: 2025年2月
**版本**: 1.0
