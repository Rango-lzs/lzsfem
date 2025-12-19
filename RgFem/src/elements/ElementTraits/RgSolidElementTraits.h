/*********************************************************************
 * \file   FEElementTraits.h
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once

#include "datastructure/Matrix.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "elements/ElementTraits/RgElementTraits.h"
#include "elements/RgElemTypeDefine.h"
#include "femcore/fem_export.h"
#include "elements/RgGaussPoint.h"

#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;
class FESolidElementShape;
class FESurfaceElementShape;


//=============================================================================
//      S O L I D   E L E M E N T
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
//! This class defines the specific traits for solid elements and serves as
//! a base class for specific solid element formulations
//
class FEM_EXPORT RgSolidElementTraits : public RgElementTraits
{
public:
    //! constructor
    RgSolidElementTraits(int ni, int ne, ElementShape es, ElementType et);

    //! initialize element traits data
    void init() override;

    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t);

    int shapeSize(int order) override;

public:
    // gauss-points
    std::vector<RgGaussPoint> gaussPoints;

    // element shape class
    FESolidElementShape* m_shape;
    // local derivatives of shape functions at gauss points
    Matrix m_Gr, m_Gs, m_Gt;
    // local second derivatives of shape functions at gauss points
    Matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};