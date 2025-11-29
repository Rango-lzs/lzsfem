# Beam 单元派生体系梳理与标准化

## 当前状态分析

### Solid 元素的派生体系（标准模型）

```
RgElement (基础类)
  └── RgSolidElement (固体元素抽象基类)
      ├── RgSolid2dElement (2D固体类)
      │   ├── RgLinearSolid2dElement (线性2D - 完全实现)
      │   │   ├── RgTri3Element
      │   │   ├── RgQuad4Element
      │   │   ├── RgQuad8Element
      │   │   └── RgTri6Element
      │   └── RgNLSolid2dElement (非线性2D - 抽象基类)
      │       ├── RgNLTri3Element
      │       └── RgNLQuad4Element
      └── RgSolid3dElement (3D固体类)
          ├── RgLinearSolid3dElement (线性3D)
          │   ├── RgTet4Element
          │   ├── RgTet10Element
          │   ├── RgHex8Element
          │   ├── RgHex20Element
          │   ├── RgWedge6Element
          │   └── RgPyramid5Element
          └── RgNLSolid3dElement (非线性3D)
              └── RgHex8GeomNLElement
```

**关键特点：**
1. ✅ 中间维度层（2D vs 3D）
2. ✅ 中间线性/非线性层
3. ✅ 具体元素类（最末层）
4. ✅ 清晰的分类逻辑

---

## Beam 元素当前派生体系（需要调整）

```
RgStructureElement (结构元素基类)
  └── RgBeamElement (梁元素抽象基类)
      ├── RgLinearBeamElement (线性梁基类)  ← 中间层
      │   ├── RgBeam2dElement (2D线性梁)
      │   └── RgBeam3dElement (3D线性梁)
      └── RgNLBeamElement (非线性梁基类)     ← 中间层
          └── RgBeam3dGeomNLElement (3D非线性梁)
```

**问题：**
1. ❌ 缺少维度层（2D vs 3D）的中间基类
2. ❌ Beam 是 1D 元素，应该有更清晰的维度划分
3. ❌ 与 Solid 体系的对称性不够好

---

## 推荐改进方案

### 方案 A: 引入 1D 维度层（最佳）

```
RgStructureElement
  └── RgBeamElement (1D梁元素基类)
      ├── RgLinearBeamElement (线性梁基类)
      │   ├── RgLinearBeam1dElement (1D线性梁基类)
      │   │   ├── RgBeam2dElement (2D/平面梁)
      │   │   └── RgBeam3dElement (3D/空间梁)
      │   └── [预留] RgLinearBeam2dElement (2D特定)
      └── RgNLBeamElement (非线性梁基类)
          ├── RgNLBeam1dElement (1D非线性梁基类)
          │   ├── RgBeam2dGeomNLElement (2D/平面非线性梁)
          │   └── RgBeam3dGeomNLElement (3D/空间非线性梁)
          └── [预留] RgNLBeam2dElement (2D特定)
```

**优点：**
- 与 Solid 的 RgSolid2dElement/RgSolid3dElement 对应
- 为未来扩展保留空间
- 逻辑更清晰

### 方案 B: 直接维度分类（折衷方案）

```
RgStructureElement
  └── RgBeamElement (1D梁元素基类)
      ├── RgLinearBeamElement
      │   ├── RgBeam2dElement (2D线性梁)
      │   └── RgBeam3dElement (3D线性梁)
      └── RgNLBeamElement
          ├── RgBeam2dGeomNLElement (2D非线性梁) ← 待创建
          └── RgBeam3dGeomNLElement (3D非线性梁) ← 已有
```

**优点：**
- 改动最小
- 具体元素类直接标明维度和线性/非线性特性
- 易于理解

---

## Shell 元素派生体系（参考）

```
RgElement
  └── RgLinearShellElement (线性壳基类)
      ├── RgShell3Element (3节点线性三角壳)
      └── RgShell4Element (4节点双线性四边形壳)
  └── RgNLShellElement (非线性壳基类)
      ├── RgShell3GeomNLElement (3节点非线性三角壳) ← 待创建
      └── RgShell4GeomNLElement (4节点非线性四边形壳) ← 待创建
```

**问题：**
- Shell 没有维度层（但壳单元本身就是2D+厚度，所以可以接受）

---

## 建议实施方案（方案 B - 改动最小）

### 需要添加的文件：

1. **RgBeam2dGeomNLElement** (2D 非线性梁)
   - 继承：RgNLBeamElement
   - DOF：4 个/节点（ux, uy, rz, theta_z）
   - 特性：平面 Timoshenko 梁，几何非线性

### 需要重命名/更新的文件：

**RgBeam3dGeomNLElement** - 现有，已继承正确的 RgNLBeamElement ✅

### 当前继承关系检查清单：

```
RgBeam2dElement
├── 继承: RgLinearBeamElement ✅
├── 维度: 2D ✅
├── 线性/非线性: Linear ✅
└── 状态: 完成

RgBeam3dElement
├── 继承: RgLinearBeamElement ✅
├── 维度: 3D ✅
├── 线性/非线性: Linear ✅
└── 状态: 完成

RgBeam3dGeomNLElement
├── 继承: RgNLBeamElement ✅
├── 维度: 3D ✅
├── 线性/非线性: NonLinear ✅
└── 状态: 完成

RgBeam2dGeomNLElement
├── 继承: RgNLBeamElement ⏳
├── 维度: 2D ⏳
├── 线性/非线性: NonLinear ⏳
└── 状态: 需要创建
```

---

## 完整的规范化派生体系（目标）

### Beam 元素完整体系

```
层级 1: RgBeamElement (基础)
        ├─ 线性/非线性分支 (层级 2)
        │  ├── RgLinearBeamElement
        │  │   ├── RgBeam2dElement (平面梁，线性)
        │  │   └── RgBeam3dElement (空间梁，线性)
        │  └── RgNLBeamElement
        │      ├── RgBeam2dGeomNLElement (平面梁，非线性)
        │      └── RgBeam3dGeomNLElement (空间梁，非线性)
        └─ [可选] 维度分支（高级）
           ├── RgLinearBeam1dElement
           └── RgNLBeam1dElement
```

### Solid 元素完整体系（参考）

```
层级 1: RgSolidElement (基础)
        ├─ 维度分支 (层级 2)
        │  ├── RgSolid2dElement
        │  │   ├── RgLinearSolid2dElement
        │  │   │   ├── RgTri3Element
        │  │   │   └── ...
        │  │   └── RgNLSolid2dElement
        │  │       ├── RgNLTri3Element
        │  │       └── ...
        │  └── RgSolid3dElement
        │      ├── RgLinearSolid3dElement
        │      │   ├── RgTet4Element
        │      │   └── ...
        │      └── RgNLSolid3dElement
        │          └── RgHex8GeomNLElement
        └─ 线性/非线性分支（已通过维度实现）
```

---

## 实施建议

### 立即行动：
1. ✅ 已完成：RgLinearBeamElement, RgNLBeamElement 创建
2. ✅ 已完成：RgBeam2dElement, RgBeam3dElement 改为继承 RgLinearBeamElement
3. ✅ 已完成：RgBeam3dGeomNLElement 改为继承 RgNLBeamElement

### 后续完善：
1. 📝 创建 RgBeam2dGeomNLElement（2D 非线性梁）
2. 📝 创建 RgShell3GeomNLElement（3D 非线性三角壳）
3. 📝 创建 RgShell4GeomNLElement（4D 非线性四边形壳）
4. 📝 [可选] 创建 RgLinearBeam1dElement 和 RgNLBeam1dElement 进一步规范化

### 命名规则总结：

**Beam 元素命名规范：**
```
Rg + [维度(Beam2d/Beam3d)] + [线性性(有=Linear/GeomNL)] + Element

示例：
- RgBeam2dElement           (2D 线性梁)
- RgBeam3dElement           (3D 线性梁)
- RgBeam2dGeomNLElement     (2D 几何非线性梁)
- RgBeam3dGeomNLElement     (3D 几何非线性梁)
```

**Solid 元素命名规范：**
```
Rg + [元素类型(Tri/Quad/Tet/Hex/etc)] + [节点数] + [线性性(无=Linear, GeomNL)] + Element

示例：
- RgTri3Element             (三角形 3节点 线性)
- RgNLTri3Element           (三角形 3节点 非线性)
- RgHex8GeomNLElement       (六面体 8节点 几何非线性)
```

---

**当前状态**: Beam 体系已 70% 规范化，Shell 体系已 50% 规范化

**建议优先级**:
1. 创建 RgBeam2dGeomNLElement (完成 Beam 体系)
2. 创建 RgShell3/4GeomNLElement (完成 Shell 体系)
3. [可选] 添加更多二次梁/壳单元
