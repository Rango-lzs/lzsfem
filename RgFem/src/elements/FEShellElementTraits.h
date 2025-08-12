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

//=============================================================================
//          T R U S S    E L E M E N T S
//
// This section defines truss elements for 3D analysis
//=============================================================================

//=============================================================================
class FETrussElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FETrussElementTraits() : FEElementTraits(NINT, NELN, FE_ELEM_TRUSS, ET_TRUSS2, FE_TRUSS) { init(); }

	void init();
};

//=============================================================================
//          D I S C R E T E    E L E M E N T S
//
// This section defines discrete elements for 3D analysis
//=============================================================================

//=============================================================================
class FEDiscreteElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FEDiscreteElementTraits() : FEElementTraits(NINT, NELN, FE_ELEM_DISCRETE, ET_DISCRETE, FE_DISCRETE) { init(); }

	void init() {}
};

//=============================================================================
//                      2 D   E L E M E N T S
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
// This class defines the traits for 2D elements and serves as a
// base class for the specific 2D element formulations.
class FE2DElementTraits: public FEElementTraits
{
public:
	FE2DElementTraits(int ni, int ne, ElementShape es, ElementType et);

	// initialization
	void init();

	// shape functions at (r,s)
	virtual void shape(double* H, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv(double* Gr, double* Gs, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) = 0;

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	Matrix Gr, Gs;
    
    // local second derivatives of shape functions at gauss points
    Matrix Grr, Gsr, Grs, Gss;
};

//=============================================================================
//
//   FE2DTri3
//   
//=============================================================================

//=============================================================================
//! Base class for linear triangles
class FE2DTri3_ : public FE2DElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	FE2DTri3_(int ni, ElementType et) : FE2DElementTraits(ni, NELN, ET_TRI3, et){}

	//! shape function at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class FE2DTri3G1 : public FE2DTri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FE2DTri3G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FE2DTri6
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles
class FE2DTri6_ : public FE2DElementTraits
{
public:
	enum { NELN = 6 };

public:
	FE2DTri6_(int ni, ElementType et) : FE2DElementTraits(ni, NELN, ET_TRI6, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//  6-node triangular element with 3-point gaussian quadrature
//
class FE2DTri6G3 : public FE2DTri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FE2DTri6G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FE2DQuad4
//   
//=============================================================================

//=============================================================================
// Base class for 4-node bilinear quadrilaterals
//
class FE2DQuad4_ : public FE2DElementTraits
{
public:
	enum { NELN = 4 };

public:
	//! constructor
	FE2DQuad4_(int ni, ElementType et) : FE2DElementTraits(ni, NELN, ET_QUAD4, et){}

	//! shape functions at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class FE2DQuad4G4 : public FE2DQuad4_
{
public:
	enum { NINT = 4 };

public:
	FE2DQuad4G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   FE2DQuad8
//   
//=============================================================================

//=============================================================================
//! Base class for 8-node quadratic quadrilaterals
//
class FE2DQuad8_ : public FE2DElementTraits
{
public:
	enum { NELN = 8 };

public:
	FE2DQuad8_(int ni, ElementType et) : FE2DElementTraits(ni, NELN, ET_QUAD8, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 9 integration points
//
class FE2DQuad8G9 : public FE2DQuad8_
{
public:
	enum { NINT = 9 };

	// constructor
	FE2DQuad8G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//
//   FE2DQuad9
//   
//=============================================================================

//=============================================================================
//! Base class for 9-node quadratic quadrilaterals
//
class FE2DQuad9_ : public FE2DElementTraits
{
public:
	enum { NELN = 9 };

public:
	FE2DQuad9_(int ni, ElementType et) : FE2DElementTraits(ni, NELN, ET_QUAD9, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 integration points
//
class FE2DQuad9G9 : public FE2DQuad9_
{
public:
	enum { NINT = 9 };

	// constructor
	FE2DQuad9G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//                      L I N E   E L E M E N T S
//
// This section defines a set of element formulations used to describe edges. 
//=============================================================================

//=============================================================================
class FELineElementTraits : public FEElementTraits
{
public:
	FELineElementTraits(int ni, int ne, ElementShape es, ElementType et);

	// initialization
	void init();

	// shape functions at r
	virtual void shape(double* H, double r) = 0;

	// shape function derivatives at (r)
	virtual void shape_deriv(double* Gr, double r) = 0;

	// shape function second derivatives at (r)
	virtual void shape_deriv2(double* Grr, double r) = 0;

public:
	std::vector<double> gr;	//!< integration point coordinates
	std::vector<double> gw;	//!< integration point weights

	// local derivatives of shape functions at gauss points
	Matrix Gr;
    
    // local second derivatives of shape functions at gauss points
    Matrix Grr;
};

//=============================================================================
//
//   FELine2_
//   
//=============================================================================

//=============================================================================
//! Base class for two-point lines
class FELine2_ : public FELineElementTraits
{
public:
	enum { NELN = 2 };

public:
	//! constructor
	FELine2_(int ni, ElementType et) : FELineElementTraits(ni, NELN, ET_LINE2, et){}

	//! shape function at (r)
	void shape(double* H, double r);

	//! shape function derivatives at (r)
	void shape_deriv(double* Gr, double r);

	//! shape function derivatives at (r)
	void shape_deriv2(double* Grr, double r);
};

//=============================================================================
class FELine2G1 : public FELine2_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FELine2G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};
