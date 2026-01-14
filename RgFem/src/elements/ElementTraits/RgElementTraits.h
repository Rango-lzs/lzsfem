/*********************************************************************
 * \file   RgElementTraits.h
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once

#include "datastructure/Matrix.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "elements/RgElemTypeDefine.h"
#include "elements/RgGaussPoint.h"
#include "femcore/fem_export.h"

#include <vector>

class NaturalCoord;

class RgElementShape;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes
//! 嵥Ԫ(Traits)͡״, ֹ򣬴洢״ڸ˹ֵ㴦ֵ
//! ElementTraitsElementShapeÿֵԪ͵ĵ
// ElementSpecify
class FEM_EXPORT RgElementTraits
    {
    public:
        //! constructor , ni ֵ ne ڵ
        RgElementTraits(int ni, int ne, ElementCategory c, ElementShape s, ElementType t);

        //! destructor
        virtual ~RgElementTraits()
        {
        }

        //! return the element class
        ElementCategory Class() const
        {
            return m_spec.eclass;
        }

        //! return the element shape
        ElementShape Shape() const
        {
            return m_spec.eshape;
        }

        //! return the element type
        ElementType Type() const
        {
            return m_spec.etype;
        }

        int shapeSize()
        {
            return m_neln;
        }

        int guassSize()
        {
            return m_nint;
        }

        RgGaussPoint gaussPoint(int i) const;

        virtual void project_to_nodes(double* ai, double* ao) const
        {
        }

        int Faces() const
        {
            return m_faces;
        }

        const Matrix& getH()
        {
            return m_H;
        }

        // values of shape functions with size N
        virtual std::vector<double> evalH(const NaturalCoord& coord) = 0;

        // values of shape function derivatives with size 3,N (2,N for 2d)
        virtual std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) = 0;

        // values of shape function second derivatives with size 6,N (3,N for 2d)
        virtual std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) = 0;

    protected:
        //! function to allocate storage for integration point data
        virtual void init() = 0;

        int m_nint;  //!< number of integration points
        int m_neln;  //!< number of element nodes

                 // gauss-points
        std::vector<RgGaussPoint> gaussPoints;

        // element shape class
        RgElementShape* m_shape = nullptr;

        Matrix m_H;              //!< shape function values at gausspoints.
                                 //!< The first index refers to the gauss-point,
                                 //!< the second index to the shape function

        FE_Element_Spec m_spec;  //!< element specs

        // number of faces of element
        int m_faces;
    };