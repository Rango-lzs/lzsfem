#pragma once

#include "elements/FEElementLibrary.h"
#include "elements/FEElementTraits.h"
#include "elements/RgElementState.h"
#include "femcore/fem_export.h"
#include "datastructure/Matrix3d.h"

#include <vector>

class FEMaterialPoint;

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
    static constexpr int MAX_NODES = 27;
    static constexpr int MAX_INTPOINTS = 27;

public:
    // Element Info
    FEElement();
    virtual ~FEElement()
    {
    }

    int getID() const;
    void setID(int n);

    int getMatID() const;
    void setMatID(int id);

    void SetLocalID(int lid)
    {
        m_loc_id = lid;
    }
    int GetLocalID() const
    {
        return m_loc_id;
    }

    virtual std::string elementType() = 0;
    virtual void setNode(FENode* n, int i);
    virtual std::vector<FENode*> getElementNodes();
    virtual FENode* giveNode(int i) const = 0;

    //! Set the type of the element and initialize the traits by type
    void SetType(int ntype)
    {
        FEElementLibrary::SetElementTraits(*this, ntype);
    }

    //! Set the traits of an element
    virtual void SetTraits(FEElementTraits* ptraits);

    //! Get the element traits
    FEElementTraits* GetTraits()
    {
        return m_pTraits;
    }

    //! return number of nodes
    int Nodes() const
    {
        return m_pTraits->m_neln;
    }

public:  //Shape function
    //! return number of integration points
    int GaussPoints() const
    {
        return m_pTraits->m_nint;
    }

    //! shape function values
    double* H(int n)
    {
        return m_pTraits->m_H[n];
    }
    const double* H(int n) const
    {
        return m_pTraits->m_H[n];
    }


public:
    //! Get the material point data
    FEMaterialPoint* GetMaterialPoint(int n)
    {
        return m_State[n];
    }

    //! set the material point data
    void SetMaterialPointData(FEMaterialPoint* pmp, int n)
    {
        pmp->m_elem = this;
        pmp->m_index = n;
        m_State[n] = pmp;
    }

    //! serialize
    //! NOTE: state data is not serialized by the element. This has to be done by the domains.
    virtual void Serialize(DumpStream& ar);

public:  // Filed evalulate
    virtual void stiffnessMatrix(Matrix3d& stiffMat) = 0;
    virtual void loadVector(std::vector<double>&) = 0;

    virtual void calcStress(FEMaterialPoint& matPt, StressTensor& stress) = 0;
    virtual void calcStrain(FEMaterialPoint& matPt, StrainTensor& strain) = 0;

    // project data to nodes, from gauss point to node
    void Extrapolation(double* ai, double* ao) const
    {
        m_pTraits->project_to_nodes(ai, ao);
    }
    void Extrapolation(vec3d<3>* ai, vec3d* ao) const
    {
        m_pTraits->project_to_nodes(ai, ao);
    }
    void Extrapolation(mat3ds* ai, mat3ds* ao) const
    {
        m_pTraits->project_to_nodes(ai, ao);
    }
    void Extrapolation(mat3d* ai, mat3d* ao) const
    {
        m_pTraits->project_to_nodes(ai, ao);
    }

    int ShapeFunctions(int order);
    double* H(int order, int n);

protected:
    int m_id;                 //!< element ID
    int m_loc_id;             //!< local ID
    int m_mat_id;             //!< material index

    std::vector<int> m_node;  //!< connectivity
    // This array stores the local node numbers, that is the node numbers
    // into the node list of a domain.
    std::vector<int> m_loc_node;  //!< local connectivity

    FEElementState m_state;
    FEElementTraits* m_pTraits;  //!< pointer to element traits
};
