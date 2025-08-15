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
//    S U R F A C E   E L E M E N T S
//
// This section defines a set of surface element formulations for use in 3D
// finite element models. For specific, these elements are used to define
// the surface of 3D volume models.
//=============================================================================

//=============================================================================
// This class defines the traits for surface elements and serves as a
// base class for the specific surface element formulations.
class FEM_EXPORT FESurfaceElementTraits : public FEElementTraits
{
public:
	FESurfaceElementTraits(int ni, int ne, ElementShape es, ElementType et);

	// initialization
	void init();

	// shape functions at (r,s)
	void shape_fnc(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);

	// shape functions at (r,s)
	void shape_fnc(int order, double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(int order, double* Gr, double* Gs, double r, double s);

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gw;

	// element shape class
	FESurfaceElementShape*				m_shape;
	std::vector<FESurfaceElementShape*>	m_shapeP; // shape classes for different order (some orders can by null)

	// local derivatives of shape functions at gauss points
	Matrix Gr, Gs;

	// local derivatives of shape functions at gauss points, for different interpolation order
	std::vector<Matrix>	Gr_p;
	std::vector<Matrix>	Gs_p;
    
    // parametric coordinates of element center
    double cr;
    double cs;
};

//=============================================================================
//
//   FEQuad4
//   
//=============================================================================

//=============================================================================
// Base class for 4-node bilinear quadrilaterals
//
class FEQuad4_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 4 };

	void init() override;

public:
	//! constructor
	FEQuad4_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_QUAD4, et){}
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class FEQuad4G4 : public FEQuad4_
{
public:
	enum { NINT = 4 };

public:
	FEQuad4G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 4-node quadrilateral elements with nodal quadrature 
class FEQuad4NI : public FEQuad4_
{
public:
	enum { NINT = 4 };

public:
	//! constructor
	FEQuad4NI();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FETri3
//   
//=============================================================================

//=============================================================================
//! Base class for linear triangles
class FETri3_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	FETri3_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_TRI3, et){}

	// initialization 
	void init() override;
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class FETri3G1 : public FETri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FETri3G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//!  3-node triangular element with 3-point gaussian quadrature
class FETri3G3 : public FETri3_
{
public:
	enum { NINT = 3 };

public:
	//! constructor
	FETri3G3();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//!  3-node triangular element with 7-point gaussian quadrature
class FETri3G7 : public FETri3_
{
public:
	enum { NINT = 7 };

public:
	//! constructor
	FETri3G7();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix	m_Ai;
};

//=============================================================================
//  3-node triangular element with nodal quadrature
class FETri3NI : public FETri3_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri3NI();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FETri6
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles
class FETri6_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 6 };

public:
	FETri6_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_TRI6, et){}

	// initialization 
	void init() override;
};

//=============================================================================
//  6-node triangular element with 3-point gaussian quadrature
//
class FETri6G3 : public FETri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri6G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 4-point gaussian quadrature
//
class FETri6G4 : public FETri6_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	FETri6G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 7-point gaussian quadrature
//
class FETri6G7 : public FETri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//  6-node triangular element with 7-point Gauss-Lobatto quadrature
//
class FETri6GL7 : public FETri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6GL7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 6-point nodal quadrature
//
class FETri6NI : public FETri6_
{
public:
	enum { NINT = 6 };

public:
	// constructor
	FETri6NI();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FETri6m
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles with modified shape functions
/*
class FETri6m_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 6 };

public:
	FETri6m_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_TRI6, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};


//=============================================================================
// 6-node triangular element (with modified shape functions)
// with 7-point gaussian quadrature
//
class FETri6mG7 : public FETri6m_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6mG7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};
*/

//=============================================================================
//
//   FETri7
//   
//=============================================================================

//=============================================================================
// Base class for 7-noded quadratic triangles
class FETri7_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 7 };

public:
	FETri7_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_TRI7, et){}

	// initialization
	void init() override;
};

//=============================================================================
//  7-node triangular element with 3-point gaussian quadrature
//
class FETri7G3 : public FETri7_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri7G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	Ai;
};

//=============================================================================
//  7-node triangular element with 4-point gaussian quadrature
//
class FETri7G4 : public FETri7_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	FETri7G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	Ai;
};


//=============================================================================
//  7-node triangular element with 7-point gaussian quadrature
//
class FETri7G7 : public FETri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri7G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//  7-node triangular element with 7-point Gauss-Lobatto quadrature
//
class FETri7GL7 : public FETri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri7GL7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	Ai;
};

//=============================================================================
//
//   FETri10
//   
//=============================================================================

//=============================================================================
// Base class for 10-noded cubic triangles
class FETri10_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 10 };

public:
	FETri10_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_TRI10, et){}

	// initialization
	void init() override;
};

//=============================================================================
//  10-node triangular element with 7-point gaussian quadrature
//
class FETri10G7 : public FETri10_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri10G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};


//=============================================================================
//  10-node triangular element with 12-point gaussian quadrature
//
class FETri10G12 : public FETri10_
{
public:
	enum { NINT = 12 };

public:
	// constructor
	FETri10G12();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//
//   FEQuad8
//   
//=============================================================================

//=============================================================================
//! Base class for 8-node quadratic quadrilaterals
//
class FEQuad8_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 8 };

public:
	FEQuad8_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_QUAD8, et) {}

	// initialization
	void init() override;
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 9 integration points
//
class FEQuad8G9 : public FEQuad8_
{
public:
	enum { NINT = 9 };

	// constructor
	FEQuad8G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 8 nodal integration points
//
class FEQuad8NI : public FEQuad8_
{
public:
    enum { NINT = 8 };
    
    // constructor
    FEQuad8NI();
    
    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;
    
private:
    Matrix	m_Ai;
};

//=============================================================================
//
//   FEQuad9
//   
//=============================================================================

//=============================================================================
//! Base class for 9-node quadratic quadrilaterals
//
class FEQuad9_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 9 };

public:
	FEQuad9_(int ni, ElementType et) : FESurfaceElementTraits(ni, NELN, ET_QUAD9, et) {}

	// initialization
	void init() override;
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 integration points
//
class FEQuad9G9 : public FEQuad9_
{
public:
	enum { NINT = 9 };

	// constructor
	FEQuad9G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 nodal integration points
//
class FEQuad9NI : public FEQuad9_
{
public:
    enum { NINT = 9 };
    
    // constructor
    FEQuad9NI();
    
    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;
    
private:
    Matrix	m_Ai;
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
