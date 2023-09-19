/*****************************************************************//**
 * \file   Element.h
 * \brief  The abstract base class of element
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef ELEMENT_H
#define ELEMENT_H

#include "Matrix.h"
#include <string>

//后续改为Parameter
//@name Input fields for general element.
//@{
#define _IFT_Element_mat "mat"
#define _IFT_Element_crosssect "crosssect"
#define _IFT_Element_nodes "nodes"
#define _IFT_Element_bodyload "bodyloads"
#define _IFT_Element_boundaryload "boundaryloads"
#define _IFT_Element_lcs "lcs"
#define _IFT_Element_partitions "partitions"
#define _IFT_Element_remote "remote"
#define _IFT_Element_activityTimeFunction "activityltf"
#define _IFT_Element_nip "nip"
//@}

class Node;
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

	//Setters and getters
	virtual std::string elementType() = 0;

	virtual void setNode(Node* n, int i);
	virtual Vector<Node*> getElementNodes();

	virtual void stiffnessMatrix(Matrix<double>&) = 0;

	virtual void loadVector(Vector<double>&) = 0;

private:

};
#endif