#pragma once

#include "datastructure/Matrix3d.h"
#include "elements/FEElementLibrary.h"
#include "elements/FEElementTraits.h"
#include "elements/FEElementState.h"
#include "femcore/fem_export.h"

#include <vector>

using NodeId = int;
using ElemId = int;
using MatId = int;

class FEMaterialPoint;
class FENode;
class FEMeshPartition;
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

class FEM_EXPORT FEElement
{
public:
    static constexpr int MAX_NODES = 28;
    static constexpr int MAX_INTPOINTS = 4;

    FEElement();
    virtual ~FEElement();

    FEElement(const FEElement& other);
    FEElement& operator=(const FEElement& other);

    // 单元编号
    int getId() const;
    void setId(int n);

    // 材料编号
    int getMatId() const;
    void setMatId(int id);

    // Get the mesh partition that contains this element
    FEMeshPartition* GetMeshPartition() const
    {
        return m_part;
    }

    // Set the mesh partition that contains this element
    void SetMeshPartition(FEMeshPartition* part)
    {
        m_part = part;
    }

    void setLocalId(int lid);
    int getLocalId() const;
    const std::vector<NodeId>& getNodeIds();

    virtual ElementType elementType() const = 0;
    
    virtual NodeId getNodeId(int idx) const;
    virtual NodeId getLocNodeId(int idx) const;

    virtual void setNode(FENode* n, int i);
    virtual FENode* getNode(int idx) const;

    //! Set the type of the element and initialize the traits by type
    void setType(int ntype)
    {
        FEElementLibrary::SetElementTraits(*this, ntype);
    }

    int getType() const;

    //Set the traits of an element
    virtual void SetTraits(FEElementTraits* ptraits);

    //Get the element traits
    FEElementTraits* GetTraits()
    {
        return m_pTraits;
    }

    //return number of nodes
    int NodeSize() const
    {
        return m_pTraits->m_neln;
    }

public: 

    //! return the element shape
    ElementShape Shape() const
    {
        return m_pTraits->Shape();
    }

    //return number of integration points
    int GaussPointSize() const;

    int ShapeFunctions(int order) const;

    // shape function values of derivations at gauss n
    double* H(int order, int n);

    //shape function values at gauss n
    //return value is array contain all the shape functions
    double* H(int n);
    const double* H(int n) const;



public:
    //Get the material point data
    FEMaterialPoint* GetMaterialPoint(int n);
    void SetMaterialPointData(FEMaterialPoint* pmp, int n);

    //! serialize
    //! NOTE: state data is not serialized by the element. This has to be done by the domains.
    virtual void Serialize(DumpStream& ar);

public:  // Filed evalulate
    virtual void stiffnessMatrix(Matrix3d& stiffMat) = 0;
    virtual void loadVector(std::vector<double>&) = 0;

    virtual void calcStress(FEMaterialPoint& matPt, StressTensor& stress) = 0;
    virtual void calcStrain(FEMaterialPoint& matPt, StrainTensor& strain) = 0;

  
protected:
    //下面的local是指在一个Domain里面的
    ElemId m_id;                     //!< element Id
    ElemId m_loc_id;                 //!< local Id
    MatId m_mat_id;                  //!< material index
    std::vector<NodeId> m_node;      //!< connectivity
    std::vector<NodeId> m_loc_node;  //!< local connectivity
    FEMeshPartition* m_part;

    FEElementState m_state;
    FEElementTraits* m_pTraits;  //!< pointer to element traits
};
