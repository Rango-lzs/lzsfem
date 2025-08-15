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

//=============================================================================
//     S H E L L   E L E M E N T S
//
// This section defines several shell formulations for use in 3D finite element
// analysis. 
//=============================================================================

//=============================================================================
// This class defines the specific for shell elements and serves as a base class
// for specific shell formulations
//
class FEShellElementTraits : public FEElementTraits
{
public:
	FEShellElementTraits(int ni, int ne, ElementShape es, ElementType et);

    void init();
    
    //! values of shape functions
    virtual void shape_fnc(double* H, double r, double s) = 0;
    
    //! values of shape function derivatives
    virtual void shape_deriv(double* Hr, double* Hs, double r, double s) = 0;
    
public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	Matrix Hr, Hs;
    
};

//=============================================================================
// 4-node quadrilateral elements
//
class FEShellQuad4_ : public FEShellElementTraits
{
public:
    enum { NELN = 4 };
    
public:
    FEShellQuad4_(int ni, ElementType et) : FEShellElementTraits(ni, NELN, ET_QUAD4, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 4-node quadrilateral elements with 4*2-point gaussian quadrature
//
class FEShellQuad4G8 : public FEShellQuad4_
{
public:
    enum { NINT = 8 };
    
public:
    FEShellQuad4G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 4-node quadrilateral elements with 4*3-point gaussian quadrature
//
class FEShellQuad4G12 : public FEShellQuad4_
{
public:
    enum { NINT = 12 };
    
public:
    FEShellQuad4G12();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements
//
class FEShellTri3_ : public FEShellElementTraits
{
public:
    enum { NELN = 3 };
    
public:
    FEShellTri3_(int ni, ElementType et) : FEShellElementTraits(ni, NELN, ET_TRI3, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 3-node triangular elements with 3*2-point gaussian quadrature
//
class FEShellTri3G6 : public FEShellTri3_
{
public:
    enum { NINT = 6 };
    
public:
    FEShellTri3G6();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements with 3*3-point gaussian quadrature
//
class FEShellTri3G9 : public FEShellTri3_
{
public:
    enum { NINT = 9 };
    
public:
    FEShellTri3G9();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements
//
class FEShellQuad8_ : public FEShellElementTraits
{
public:
    enum { NELN = 8 };
    
public:
    FEShellQuad8_(int ni, ElementType et) : FEShellElementTraits(ni, NELN, ET_QUAD8, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 8-node quadrilateral elements with 9*2-point gaussian quadrature
//
class FEShellQuad8G18 : public FEShellQuad8_
{
public:
    enum { NINT = 18 };
    
public:
    FEShellQuad8G18();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements with 9*3-point gaussian quadrature
//
class FEShellQuad8G27 : public FEShellQuad8_
{
public:
    enum { NINT = 27 };
    
public:
    FEShellQuad8G27();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements
//
class FEShellTri6_ : public FEShellElementTraits
{
public:
    enum { NELN = 6 };
    
public:
    FEShellTri6_(int ni, ElementType et) : FEShellElementTraits(ni, NELN, ET_TRI6, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 6-node triangular elements with 7*2-point gaussian quadrature
//
class FEShellTri6G14 : public FEShellTri6_
{
public:
    enum { NINT = 14 };
    
public:
    FEShellTri6G14();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements with 7*3-point gaussian quadrature
//
class FEShellTri6G21 : public FEShellTri6_
{
public:
    enum { NINT = 21 };
    
public:
    FEShellTri6G21();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

