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

class FEMaterialPoint;
class FENode;

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

    static constexpr int MAX_NODES = 28;

    FEElement();
    virtual ~FEElement();

    //��Ԫ���
    int getID() const;
    void setID(int n);

    //���ϱ��
    int getMatID() const;
    void setMatID(int id);

    void setLocalID(int lid);
    int getLocalID() const;

    virtual ElementType elementType() = 0;
    virtual void setNode(FENode* n, int i);
    virtual const std::vector<FENode*>& getNodes();
    virtual FENode* giveNode(int i) const = 0;

    //! Set the type of the element and initialize the traits by type
    void setType(int ntype)
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
    int NodeSize() const
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
        return m_state[n];
    }

    //! set the material point data
    void SetMaterialPointData(FEMaterialPoint* pmp, int n);
        /*  {
              pmp->m_elem = this;
              pmp->m_index = n;
              m_state[n] = pmp;
          }*/

    //! serialize
    //! NOTE: state data is not serialized by the element. This has to be done by the domains.
    virtual void Serialize(DumpStream& ar);

public:  // Filed evalulate
    virtual void stiffnessMatrix(Matrix3d& stiffMat) = 0;
    virtual void loadVector(std::vector<double>&) = 0;

    virtual void calcStress(FEMaterialPoint& matPt, StressTensor& stress) = 0;
    virtual void calcStrain(FEMaterialPoint& matPt, StrainTensor& strain) = 0;

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
