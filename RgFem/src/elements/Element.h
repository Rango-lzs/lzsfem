/*****************************************************************/ /**
                                                                     * \file   Element.h
                                                                     * \brief  The abstract base class of element
                                                                     *
                                                                     * \author Leizs
                                                                     * \date   September 2023
                                                                     *********************************************************************/

#ifndef ELEMENT_H
#define ELEMENT_H

#include "floatmatrix.h"

#include <string>

// 后续改为Parameter
//@name Input fields for general element.
//@{
#define _IFT_Element_mat                  "mat"
#define _IFT_Element_crosssect            "crosssect"
#define _IFT_Element_nodes                "nodes"
#define _IFT_Element_bodyload             "bodyloads"
#define _IFT_Element_boundaryload         "boundaryloads"
#define _IFT_Element_lcs                  "lcs"
#define _IFT_Element_partitions           "partitions"
#define _IFT_Element_remote               "remote"
#define _IFT_Element_activityTimeFunction "activityltf"
#define _IFT_Element_nip                  "nip"
//@}

namespace fem
{
    class Node;
    class MaterialPoint;
    /**
     * @~English
     * @brief brief-description-about-Element .
     * @
     *
     * @~Chinese
     * @brief brief-description-about-Element.
     * Tasks:
     *	单元相关的数据，节点，材料等
     *	计算单元刚度矩阵，载荷向量
     *	计算单元应力，应变
     *	结果输出
     *
     */
    class Element
    {
    public:
        Element() = default;
        ~Element() = default;

        // Setters and getters
        virtual std::string elementType() = 0;

        virtual void setNode(Node* n, int i);
        virtual std::vector<Node*> getElementNodes();
        virtual Node* giveNode(int i) const = 0;

        virtual void stiffnessMatrix(FloatMatrix& stiffMat) = 0;
        virtual void loadVector(std::vector<double>&) = 0;

        virtual void calcStress(MaterialPoint& matPt, FloatMatrix& stress) = 0;
        virtual void calcStrain(MaterialPoint& matPt, FloatMatrix& strain) = 0;

    private:
    };
}  // namespace fem
#endif

