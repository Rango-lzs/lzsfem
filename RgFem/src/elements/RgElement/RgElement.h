#pragma once

#include "datastructure/Matrix3d.h"
#include "elements/RgElementLibrary.h"
#include "elements/ElementTraits/RgElementTraits.h"
#include "elements/RgElementState.h"
#include "femcore/fem_export.h"

#include <vector>
#include <Eigen/Dense>

using NodeId = int;
using ElemId = int;
using MatId = int;

// Type aliases for matrices and vectors
using StressTensor = Matrix3ds;
using StrainTensor = Matrix3ds;

class FEMaterialPoint;
class FENode;
class RgDomain;
class DumpStream;

/**
 *@~English
 * @brief brief - description - about - Element .
 * @
 *
 *@~Chinese
 * @brief brief - description - about - Element.
 * Tasks:
 * 单元相关的数据，节点，材料等  ： 属性数据
 * 计算单元刚度矩阵，载荷向量	： 物理特性
 * 计算单元应力，应变等结果		： 结果数据
 * 结果输出
 */

/* 单元如何分类
 * 1、实体单元(连续体单元)，3D Solid，2D Plane,
 * 2、结构单元, Shell, Beam, Truss
 * 3、连接单元, Spring, Cohesive
 * 4、
 */

class FEM_EXPORT RgElement
{
public:
    static constexpr int MAX_NODES = 28;
    static constexpr int MAX_INTPOINTS = 4;

    RgElement();
    virtual ~RgElement(){};
    RgElement(const RgElement& other) = default;
    RgElement& operator=(const RgElement& other) = default;

    // --- Element Identification ---
    int getId() const;
    void setId(int n);

    int getMatId() const;
    void setMatId(int id);

    void setLocalId(int lid);
    int getLocalId() const;

    // --- Mesh Partition ---
    RgDomain* getDomain() const;
    void setDomain(RgDomain* dom);

    // --- Node Connectivity ---
    const std::vector<NodeId>& getNodeIds() const;
    virtual NodeId getNodeId(int idx) const;
    virtual NodeId getLocNodeId(int idx) const;
    void setNodeId(int idx, int id);

    virtual void setNode(FENode* n, int i);
    virtual FENode* getNode(int idx) const;

    // --- Element Properties ---
    virtual ElementType elementType() const;
    virtual ElementCategory Class() const;
    virtual ElementShape Shape() const;

    // --- Traits Management ---
    virtual void initTraits();
    RgElementTraits* getTraits();

    int NodeSize() const;
    int GaussPointSize() const;
    int ShapeFunctions(int order) const;

    // --- Material Point Data ---
    RgMaterialPoint* getMaterialPoint(int n);
    void setMaterialPointData(RgMaterialPoint* pmp, int n);

    // --- Serialization ---
    virtual void Serialize(DumpStream& ar);

    // --- Field Evaluation ---
    Vector3d Evaluate(Vector3d* value, int iGauss);

    // --- Core Finite Element Methods ---
    virtual void calculateStiffnessMatrix(Matrix& K) const = 0;
    virtual void calculateMassMatrix(Matrix& M) const = 0;
    virtual void calculateDampingMatrix(Matrix& C) const = 0;
    virtual void calculateInternalForceVector(std::vector<double>& F) const = 0;

    virtual void calculateStress(FEMaterialPoint& matPt, StressTensor& stress) = 0;
    virtual void calculateStrain(FEMaterialPoint& matPt, StrainTensor& strain) = 0;

    bool isActive() const;
    void ClearData();

protected:
    // Element data
    ElemId m_id;                     //!< element Id
    ElemId m_loc_id;                 //!< local Id in the domain
    MatId m_mat_id;                  //!< material index  
    RgDomain* m_part;

    RgElementState m_state;
    RgElementTraits* m_pTraits;      //!< pointer to element traits
    
    std::vector<NodeId> m_node;      //!< connectivity
    std::vector<NodeId> m_loc_node;  //!< local connectivity
};

