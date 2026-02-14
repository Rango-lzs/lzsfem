# RgLoad 载荷系统使用指南

## 目录
1. [概述](#概述)
2. [载荷类型](#载荷类型)
3. [Abaqus格式解析](#abaqus格式解析)
4. [使用示例](#使用示例)
5. [API参考](#api参考)

---

## 概述

本载荷系统提供完整的有限元载荷管理功能，包括：

- ✅ **节点力载荷** (Concentrated loads - *Cload)
- ✅ **表面载荷** (Pressure/Traction - *Dload P/TRVEC)
- ✅ **体力载荷** (Gravity/Centrifugal - *Dload GRAV/CENTRIF)
- ✅ **力矩载荷** (Moments)
- ✅ **Abaqus输入文件解析**
- ✅ **时间变化载荷**（通过载荷控制器）

---

## 载荷类型

### 类层次结构

```
FEObjectBase
    └── RgLoad (基类)
            ├── RgNodalLoad (节点力)
            ├── RgSurfaceLoad (表面载荷)
            ├── RgBodyLoad (体力)
            └── RgMomentLoad (力矩)
```

### 1. RgNodalLoad - 节点力载荷

**用途**: 集中力，施加在节点上

**属性**:
- `m_nodeSet`: 施加载荷的节点集
- `m_force`: 力向量 (Vector3d)
- `m_dof`: 单个DOF索引 (-1表示向量载荷)

**使用场景**:
- 点载荷
- 集中力
- 边界反力检查

**示例**:
```cpp
// 单个DOF力
RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
    fem, nodeSet, 0, 100.0);  // 100N in X direction

// 向量力
Vector3d F(100, -50, 0);
RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(fem, nodeSet, F);
```

### 2. RgSurfaceLoad - 表面载荷

**用途**: 分布载荷，施加在表面上

**类型**:
- `PRESSURE` - 法向压力
- `TRACTION` - 切向牵引力
- `GENERAL_TRACTION` - 一般牵引力（有分量）

**属性**:
- `m_surface/m_facetSet`: 施加载荷的表面
- `m_loadType`: 载荷类型
- `m_traction`: 牵引力向量
- `m_bFollower`: 是否为跟随载荷

**跟随载荷**:
- `true`: 载荷方向随表面变形更新（如内压）
- `false`: 载荷方向固定（死载荷）

**示例**:
```cpp
// 压力载荷
RgSurfaceLoad* load = RgLoadFactory::CreatePressure(
    fem, facetSet, 1.0e6);  // 1 MPa

// 牵引力
Vector3d traction(100, 0, 0);
RgSurfaceLoad* load = RgLoadFactory::CreateTraction(
    fem, surface, traction);

// 跟随压力（气球膨胀）
load->SetFollower(true);
```

### 3. RgBodyLoad - 体力载荷

**用途**: 施加在整个体积上的力

**类型**:
- `GRAVITY` - 重力
- `CENTRIFUGAL` - 离心力
- `CONSTANT` - 恒定体力

**属性**:
- `m_bodyLoadType`: 体力类型
- `m_force`: 力向量（重力）
- `m_axis`: 旋转轴（离心）
- `m_origin`: 旋转中心（离心）
- `m_omega`: 角速度（离心）

**示例**:
```cpp
// 重力
Vector3d g(0, 0, -9.81);
RgBodyLoad* load = RgLoadFactory::CreateGravity(fem, g);

// 离心力
Vector3d axis(0, 0, 1);  // Z轴
Vector3d origin(0, 0, 0);
double omega = 100.0;     // 100 rad/s
RgBodyLoad* load = RgLoadFactory::CreateCentrifugal(
    fem, axis, origin, omega);
```

### 4. RgMomentLoad - 力矩载荷

**用途**: 集中力矩

**属性**:
- `m_nodeSet`: 施加力矩的节点集
- `m_moment`: 力矩向量

**示例**:
```cpp
Vector3d M(0, 0, 100);  // 100 N⋅m about Z
RgMomentLoad* load = RgLoadFactory::CreateMoment(fem, nodeSet, M);
```

---

## Abaqus格式解析

### *Cload - 集中载荷

#### 格式
```
** Name: Load-1 Type: Concentrated force
*Cload
NodeSetName, DOF, Magnitude
```

#### DOF编号 (Abaqus约定)
- 1, 2, 3 = X, Y, Z 方向的力
- 4, 5, 6 = X, Y, Z 方向的力矩

#### 示例

**例1**: 单向力
```
*Cload
LoadPoint, 1, 100.0    # 100N in X direction
```

**例2**: 多向力
```
*Cload
LoadPoint, 1, 100.0    # X force
LoadPoint, 2, -50.0    # Y force
LoadPoint, 3, 75.0     # Z force
```

**例3**: 力矩
```
*Cload
MomentPoint, 4, 50.0   # Moment about X
```

### *Dload - 分布载荷

#### 格式1: 压力
```
*Dload
SurfaceName, P, Magnitude
```

**示例**:
```
*Dload
TopFace, P, 1.0e6      # 1 MPa pressure
```

#### 格式2: 牵引力
```
*Dload
SurfaceName, TRVEC, Magnitude, x, y, z
```

**示例**:
```
*Dload
SideFace, TRVEC, 1000.0, 1.0, 0.0, 0.0  # 1000 Pa in X direction
```

#### 格式3: 重力
```
*Dload
ElementSetName, GRAV, Magnitude, x, y, z
```

**参数**:
- `Magnitude`: 重力加速度大小 (m/s²)
- `x, y, z`: 方向（会被归一化）

**示例**:
```
*Dload
AllElements, GRAV, 9.81, 0.0, 0.0, -1.0  # g in -Z direction
```

#### 格式4: 离心载荷
```
*Dload
ElementSetName, CENTRIF, Omega^2, x0, y0, z0, x1, y1, z1
```

**参数**:
- `Omega^2`: 角速度的平方
- `x0, y0, z0`: 旋转轴上的点（原点）
- `x1, y1, z1`: 旋转轴方向

**示例**:
```
*Dload
AllElements, CENTRIF, 10000.0, 0, 0, 0, 0, 0, 1
# omega = sqrt(10000) = 100 rad/s, 绕Z轴旋转
```

### 使用解析器

```cpp
// 1. 创建解析器
AbaqusLoadParser parser(fem);

// 2. 解析输入
for (size_t i = 0; i < lines.size(); ++i)
{
    if (lines[i].find("*Cload") == 0)
    {
        int consumed = parser.ParseCload(lines, i);
        i += consumed - 1;
    }
    else if (lines[i].find("*Dload") == 0)
    {
        int consumed = parser.ParseDload(lines, i);
        i += consumed - 1;
    }
}

// 3. 创建载荷
parser.CreateLoads();

// 4. 查看解析结果
const auto& loadData = parser.GetLoadData();
for (const auto& data : loadData)
{
    std::cout << "Load: " << data.name 
              << ", Type: " << data.loadType << "\n";
}
```

---

## 使用示例

### 示例1: 悬臂梁端部载荷

```cpp
void CantileverBeam(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 端部向下的力
    FENodeSet* tipNodes = mesh.FindNodeSet("BeamTip");
    RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
        fem, tipNodes, 2, -1000.0);  // -1000N in Z
    load->SetName("TipLoad");
    fem->AddModelLoad(load);
}
```

### 示例2: 压力容器

```cpp
void PressureVessel(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 内压
    FEFacetSet* innerSurface = mesh.FindFacetSet("InnerSurface");
    RgSurfaceLoad* pressure = RgLoadFactory::CreatePressure(
        fem, innerSurface, 5.0e6);  // 5 MPa
    pressure->SetFollower(true);     // 跟随载荷
    pressure->SetName("InternalPressure");
    fem->AddModelLoad(pressure);
}
```

### 示例3: 旋转机械

```cpp
void RotatingDisk(FEModel* fem)
{
    // 离心载荷
    Vector3d axis(0, 0, 1);    // Z轴
    Vector3d origin(0, 0, 0);
    double rpm = 3000.0;
    double omega = rpm * 2.0 * M_PI / 60.0;  // Convert to rad/s
    
    RgBodyLoad* load = RgLoadFactory::CreateCentrifugal(
        fem, axis, origin, omega);
    load->SetName("CentrifugalLoad");
    fem->AddModelLoad(load);
}
```

### 示例4: 时变载荷

```cpp
void TimeVaryingLoad(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    FENodeSet* nodeSet = mesh.FindNodeSet("LoadPoint");
    
    // 创建载荷曲线
    RgLoadCurve* lc = new RgLoadCurve();
    lc->AddPoint(0.0, 0.0);     // 初始
    lc->AddPoint(1.0, 1.0);     // 加载
    lc->AddPoint(2.0, 1.0);     // 保持
    lc->AddPoint(3.0, 0.0);     // 卸载
    lc->SetInterpolation(RgLoadController::INTERP_SMOOTH);
    fem->AddLoadController(lc);
    
    // 创建载荷
    RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
        fem, nodeSet, 2, 5000.0);  // 5000N max
    load->SetLoadController(lc);
    load->SetName("RampedLoad");
    fem->AddModelLoad(load);
}
```

### 示例5: 组合载荷

```cpp
void CombinedLoading(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 1. 重力（恒定）
    Vector3d g(0, 0, -9.81);
    RgBodyLoad* gravity = RgLoadFactory::CreateGravity(fem, g);
    gravity->SetName("Gravity");
    fem->AddModelLoad(gravity);
    
    // 2. 点载荷（时变）
    FENodeSet* loadPoint = mesh.FindNodeSet("LoadPoint");
    RgLoadCurve* lc = new RgLoadCurve();
    lc->AddPoint(0.0, 0.0);
    lc->AddPoint(1.0, 1.0);
    fem->AddLoadController(lc);
    
    RgNodalLoad* pointLoad = RgLoadFactory::CreateNodalLoad(
        fem, loadPoint, 1, 1000.0);
    pointLoad->SetLoadController(lc);
    pointLoad->SetName("PointLoad");
    fem->AddModelLoad(pointLoad);
    
    // 3. 压力（恒定）
    FEFacetSet* topFace = mesh.FindFacetSet("TopFace");
    RgSurfaceLoad* pressure = RgLoadFactory::CreatePressure(
        fem, topFace, 0.1e6);
    pressure->SetName("Pressure");
    fem->AddModelLoad(pressure);
}
```

---

## API参考

### RgLoadFactory - 工厂类

```cpp
// 创建节点力（单DOF）
static RgNodalLoad* CreateNodalLoad(
    FEModel* fem, FENodeSet* nodeSet, int dof, double magnitude);

// 创建节点力（向量）
static RgNodalLoad* CreateNodalLoad(
    FEModel* fem, FENodeSet* nodeSet, const Vector3d& force);

// 创建压力载荷
static RgSurfaceLoad* CreatePressure(
    FEModel* fem, FEFacetSet* facetSet, double pressure);

// 创建牵引力载荷
static RgSurfaceLoad* CreateTraction(
    FEModel* fem, FESurface* surface, const Vector3d& traction);

// 创建重力载荷
static RgBodyLoad* CreateGravity(
    FEModel* fem, const Vector3d& g);

// 创建离心载荷
static RgBodyLoad* CreateCentrifugal(
    FEModel* fem, const Vector3d& axis, 
    const Vector3d& origin, double omega);

// 创建力矩载荷
static RgMomentLoad* CreateMoment(
    FEModel* fem, FENodeSet* nodeSet, const Vector3d& moment);
```

### AbaqusLoadParser

```cpp
// 构造函数
AbaqusLoadParser(FEModel* fem);

// 解析 *Cload
int ParseCload(const std::vector<std::string>& lines, int startIdx);

// 解析 *Dload
int ParseDload(const std::vector<std::string>& lines, int startIdx);

// 创建载荷
bool CreateLoads();

// 获取解析数据
const std::vector<LoadData>& GetLoadData() const;
```

---

## 集成到求解器

### 初始化

```cpp
bool FEModel::Init()
{
    // ... 其他初始化 ...
    
    // 初始化载荷
    int nloads = ModelLoads();
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(ModelLoad(i));
        if (!load->Init())
        {
            feLogError("Failed to initialize load %d", i);
            return false;
        }
    }
    
    // 激活载荷
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(ModelLoad(i));
        load->Activate();
    }
    
    return true;
}
```

### 时间步更新

```cpp
void FEModel::Update()
{
    // 评估载荷控制器
    const FETimeInfo& tp = GetTime();
    EvaluateLoadControllers(tp.currentTime);
    
    // 更新载荷
    int nloads = ModelLoads();
    for (int i = 0; i < nloads; ++i)
    {
        RgLoad* load = dynamic_cast<RgLoad*>(ModelLoad(i));
        if (load->IsActive())
        {
            load->Update();
        }
    }
}
```

### 装配到整体载荷向量

节点载荷会直接设置到节点的载荷数组中：
```cpp
// 在 RgNodalLoad::Update() 中
node.set_load(dof, magnitude * loadController->Value());
```

表面载荷和体力载荷需要在单元装配时处理：
```cpp
// 在域的 InternalForces 或 BodyForce 中
void RgSolidDomain::BodyForce(FEGlobalVector& R, RgBodyLoad& bf)
{
    for (int i = 0; i < Elements(); ++i)
    {
        // 计算单元体力
        Vector3d f = bf.EvaluateForce(materialPoint);
        // 装配到整体向量
    }
}
```

---

## 常见问题

### Q1: 如何施加分布载荷到表面？

A: 使用 `RgSurfaceLoad`，可以是压力或牵引力。

```cpp
FEFacetSet* facetSet = mesh.FindFacetSet("Surface");
RgSurfaceLoad* load = RgLoadFactory::CreatePressure(fem, facetSet, 1.0e6);
```

### Q2: 跟随载荷是什么？

A: 跟随载荷的方向随着结构变形而更新。例如内压始终垂直于变形后的表面。

```cpp
load->SetFollower(true);
```

### Q3: 如何实现循环载荷？

A: 使用 `RgSineController`。

```cpp
RgSineController* lc = new RgSineController();
lc->SetParameters(amplitude, frequency, phase, offset);
load->SetLoadController(lc);
```

### Q4: 重力和体力有什么区别？

A: 重力是特殊的体力，按密度缩放。体力可以是任意的体积力。

---

## 总结

RgLoad 载荷系统提供：
- ✅ 完整的载荷类型
- ✅ Abaqus格式兼容
- ✅ 时间变化支持
- ✅ 跟随载荷
- ✅ 工厂模式
- ✅ 易于扩展

可直接集成到FEM求解器使用！
