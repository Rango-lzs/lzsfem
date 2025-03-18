#pragma once

#include "datastructure/Matrix3d.h"
#include "elements/FEElementLibrary.h"
#include "elements/FEElementTraits.h"
#include "elements/RgElementState.h"
#include "femcore/fem_export.h"

#include <vector>

using NodeId = int;
using ElemId = int;
using MatId = int;

class RgMaterialPoint;

/**
 *@~English
 * @brief brief - description - about - Element .
 * @
 *
 *@~Chinese
 * @brief brief - description - about - Element.
 * Tasks:
 * ��Ԫ��ص����ݣ��ڵ㣬���ϵ�  �� ��������
 * ���㵥Ԫ�նȾ����غ�����	�� ��������
 * ���㵥ԪӦ����Ӧ��Ƚ��		�� �������
 * ������
 */

/* ��Ԫ��η���
 * 1��ʵ�嵥Ԫ(�����嵥Ԫ)��3D Solid��2D Plane,
 * 2���ṹ��Ԫ, Shell, Beam, Truss
 * 3�����ӵ�Ԫ, Spring, Cohesive
 * 4��
 */

class FEM_EXPORT FEElement
{
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

    void SetLocalID(int lid);
    int GetLocalID() const;

    virtual ElementType elementType() = 0;
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

public:  // Shape function
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
    ElemId m_id;                 //!< element ID
    ElemId m_loc_id;             //!< local ID
    MatId m_mat_id;              //!< material index
    std::vector<NodeId> m_node;  //!< connectivity

    FEElementState m_state;
    FEElementTraits* m_pTraits;  //!< pointer to element traits
};
