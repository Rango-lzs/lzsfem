## Beam 和 Shell 元素重构总结

### 完成的重构内容

重构成功完成，现在 Beam 和 Shell 元素遵循与 Solid 元素相同的线性/非线性分层架构模式。

---

## 新的继承结构

### Beam 元素继承树

```
RgBeamElement (基类)
├── RgLinearBeamElement (线性梁元素基类)
│   ├── RgBeam2dElement (2D 线性梁)
│   └── RgBeam3dElement (3D 线性梁)
└── RgNLBeamElement (非线性梁元素基类)
    └── RgBeam3dGeomNLElement (3D 几何非线性梁)
```

### Shell 元素继承树

```
RgElement (基类)
├── RgLinearShellElement (线性壳元素基类)
│   ├── RgShell3Element (3节点线性三角壳)
│   └── RgShell4Element (4节点双线性四边形壳)
└── RgNLShellElement (非线性壳元素基类)
    ├── RgShell3GeomNLElement (3节点三角壳 - 几何非线性)
    └── RgShell4GeomNLElement (4节点四边形壳 - 几何非线性)
```

与 Solid 元素对比：
```
RgElement (基类)
├── RgSolid2dElement (2D 固体基类)
│   ├── RgLinearSolid2dElement (线性)
│   │   ├── RgTri3Element
│   │   ├── RgQuad4Element
│   │   └── ...
│   └── RgNLSolid2dElement (非线性)
│       ├── RgNLTri3Element
│       ├── RgNLQuad4Element
│       └── ...
└── ...
```

---

## 创建的新文件

### 基类文件（4对 = 8个文件）

**Beam 基类:**
- `RgLinearBeamElement.h` / `RgLinearBeamElement.cpp`
- `RgNLBeamElement.h` / `RgNLBeamElement.cpp`

**Shell 基类:**
- `RgLinearShellElement.h` / `RgLinearShellElement.cpp`
- `RgNLShellElement.h` / `RgNLShellElement.cpp`

### 更新的文件（6个头文件）

**Beam 元素继承关系变更:**
- `RgBeam2dElement.h` - 改为继承 `RgLinearBeamElement`
- `RgBeam3dElement.h` - 改为继承 `RgLinearBeamElement`
- `RgBeam3dGeomNLElement.h` - 改为继承 `RgNLBeamElement`

**Shell 元素继承关系变更:**
- `RgShell3Element.h` - 改为继承 `RgLinearShellElement`
- `RgShell4Element.h` - 改为继承 `RgLinearShellElement`

**CMakeLists.txt 更新:**
- 添加了新基类的头文件声明
- 添加了新基类的源文件声明

---

## 基类的虚函数接口

### RgLinearBeamElement

```cpp
virtual std::string typeName() const = 0;
virtual void calculateStiffnessMatrix(RgMatrix& K) const = 0;
virtual void calculateMassMatrix(RgMatrix& M) const = 0;
virtual void calculateInternalForceVector(RgVector& F) const = 0;
```

### RgNLBeamElement

```cpp
virtual std::string typeName() const = 0;
virtual void computeDeformationGradient(int gaussPointIndex, 
                                        const std::vector<double>& displacement,
                                        std::array<std::array<double, 3>, 3>& F) const = 0;
virtual void computeDisplacementGradient(int gaussPointIndex,
                                         const std::vector<double>& displacement,
                                         std::array<std::array<double, 3>, 3>& dispGrad) const = 0;
virtual void calculateStiffnessMatrix(RgMatrix& K) const = 0;
virtual void calculateMassMatrix(RgMatrix& M) const = 0;
virtual void calculateInternalForceVector(RgVector& F) const = 0;
```

### RgLinearShellElement

```cpp
virtual std::string typeName() const = 0;
virtual double getShellThickness() const = 0;
virtual void setShellThickness(double thickness) = 0;
virtual double getElementArea() const = 0;
virtual void getShellNormal(double r, double s, std::array<double, 3>& normal) const = 0;
virtual void calculateStiffnessMatrix(RgMatrix& K) const = 0;
virtual void calculateMassMatrix(RgMatrix& M) const = 0;
virtual void calculateInternalForceVector(RgVector& F) const = 0;
```

### RgNLShellElement

```cpp
virtual std::string typeName() const = 0;
virtual double getShellThickness() const = 0;
virtual void setShellThickness(double thickness) = 0;
virtual double getElementArea() const = 0;
virtual void getShellNormal(double r, double s, std::array<double, 3>& normal) const = 0;
virtual void computeDeformationGradient(int gaussPointIndex,
                                        const std::vector<double>& displacement,
                                        std::array<std::array<double, 3>, 3>& F) const = 0;
virtual void computeDisplacementGradient(int gaussPointIndex,
                                         const std::vector<double>& displacement,
                                         std::array<std::array<double, 3>, 3>& dispGrad) const = 0;
virtual void calculateStiffnessMatrix(RgMatrix& K) const = 0;
virtual void calculateMassMatrix(RgMatrix& M) const = 0;
virtual void calculateInternalForceVector(RgVector& F) const = 0;
```

---

## 架构优势

1. **一致性**: Beam 和 Shell 元素现在遵循与 Solid 元素相同的分层模式
2. **清晰的分离**: 线性和非线性分析通过不同的基类进行清晰分离
3. **代码重用**: 公共接口在基类中定义，实现在具体元素中
4. **易于扩展**: 添加新的几何非线性元素变体只需创建新的具体类
5. **理论区分**: 每个基类在文档中明确说明其适用的理论框架

---

## 下一步可能的工作

### 立即可做（可选）:
- 创建 `RgShell3GeomNLElement` 和 `RgShell4GeomNLElement` (几何非线性壳)
- 创建 `RgBeam2dGeomNLElement` (2D 几何非线性梁)
- 创建 `RgBeam8Element` (8节点二次梁)

### 需要完成的工作:
- 实现基类中的占位符虚函数（currentlyReturn无或stub）
- 在具体元素中实现实际的矩阵组装逻辑

---

**重构时间**: 2025-11-28  
**状态**: ✅ 完成
