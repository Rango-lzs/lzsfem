/*********************************************************************
 * \file   FEElementTraits.h
 * \brief  
 * 
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once

#include "femcore/fem_export.h"

#include "datastructure/Matrix.h"
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix3d.h"
#include "elements/RgElemTypeDefine.h"
#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;
class FESolidElementShape;
class FESurfaceElementShape;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes
//! ���嵥Ԫ������(Traits)��������͡���״������, ���ֹ��򣬴洢��״�����ڸ�˹���ֵ㴦��ֵ
//! ElementTraits ��ÿ����Ԫ����һ��ʵ���� ��ElementShape�����ÿ����Ԫ��ĵ���
//ElementSpecify
class FEM_EXPORT FEElementTraits
{
public:
	//! constructor , ni ���ֵ����� ne �ڵ���
	FEElementTraits(int ni, int ne, ElementCategory c, ElementShape s, ElementType t);

	//! destructor
	virtual ~FEElementTraits(){}

	//! return the element class
	ElementCategory Class() const { return m_spec.eclass; }

	//! return the element shape
	ElementShape Shape() const { return m_spec.eshape; }

	//! return the element type
	ElementType Type() const { return m_spec.etype; }

	// ���ֵ��ֵ��嵽�ڵ�
	virtual void project_to_nodes(double* ai, double* ao) const {}
	virtual void project_to_nodes(Vector3d*  ai, Vector3d*  ao) const;
	virtual void project_to_nodes(Matrix3ds* ai, Matrix3ds* ao) const;
	virtual void project_to_nodes(Matrix3d*  ai, Matrix3d*  ao) const;

	virtual int ShapeFunctions(int order) { return m_neln; }

	int Faces() const { return m_faces; }

public:
	int m_nint;	//!< number of integration points
	int	m_neln;	//!< number of element nodes

	Matrix m_H;	//!< shape function values at gausspoints.
				//!< The first index refers to the gauss-point,
				//!< the second index to the shape function

	FE_Element_Spec	m_spec;	//!< element specs

protected:
	// number of faces of element
	int	m_faces;

	//! function to allocate storage for integration point data
	virtual void init() = 0;
};
