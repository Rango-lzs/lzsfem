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
#include "elements/RgGaussPoint.h"
#include <vector>
#include "RgElementTraits.h"

class RgElement;
class RgSurfaceElementShape;

class NaturalCoord;


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
class FEM_EXPORT RgSurfaceElementTraits : public RgElementTraits
{
public:
	RgSurfaceElementTraits(int ni, int ne, ElementShape es, ElementType et);

	// initialization
	void init() override;

    // values of shape functions with size N
    virtual std::vector<double> evalH(const NaturalCoord& coord);

    // values of shape function derivatives with size 3,N (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord);

    // values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord);

public:
	// gauss-points
	std::vector<RgGaussPoint> gaussPoints;

	// local derivatives of shape functions at gauss points
	Matrix Gr, Gs;

	// local second derivatives of shape functions at gauss points
	Matrix Grr, Gsr, Gss;
};

//=============================================================================
//
//   FEQuad4
//   
//=============================================================================

//=============================================================================
// Base class for 4-node bilinear quadrilaterals
//
class RgQuad4_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 4 };

	void init() override;

public:
	//! constructor
	RgQuad4_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_QUAD4, et){}
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class RgQuad4G4 : public RgQuad4_
{
public:
	enum { NINT = 4 };

public:
	RgQuad4G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 4-node quadrilateral elements with nodal quadrature 
class RgQuad4NI : public RgQuad4_
{
public:
	enum { NINT = 4 };

public:
	//! constructor
	RgQuad4NI();

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
class RgTri3_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	RgTri3_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_TRI3, et){}

	// initialization 
	void init() override;
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class RgTri3G1 : public RgTri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	RgTri3G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//!  3-node triangular element with 3-point gaussian quadrature
class RgTri3G3 : public RgTri3_
{
public:
	enum { NINT = 3 };

public:
	//! constructor
	RgTri3G3();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//!  3-node triangular element with 7-point gaussian quadrature
class RgTri3G7 : public RgTri3_
{
public:
	enum { NINT = 7 };

public:
	//! constructor
	RgTri3G7();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix	m_Ai;
};

//=============================================================================
//  3-node triangular element with nodal quadrature
class RgTri3NI : public RgTri3_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	RgTri3NI();

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
class RgTri6_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 6 };

public:
	RgTri6_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_TRI6, et){}

	// initialization 
	void init() override;
};

//=============================================================================
//  6-node triangular element with 3-point gaussian quadrature
//
class RgTri6G3 : public RgTri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	RgTri6G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 4-point gaussian quadrature
//
class RgTri6G4 : public RgTri6_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	RgTri6G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 7-point gaussian quadrature
//
class RgTri6G7 : public RgTri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	RgTri6G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//  6-node triangular element with 7-point Gauss-Lobatto quadrature
//
class RgTri6GL7 : public RgTri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	RgTri6GL7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  6-node triangular element with 6-point nodal quadrature
//
class RgTri6NI : public RgTri6_
{
public:
	enum { NINT = 6 };

public:
	// constructor
	RgTri6NI();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// Base class for 7-noded quadratic triangles
class RgTri7_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 7 };

public:
	RgTri7_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_TRI7, et){}

	// initialization
	void init() override;
};

//=============================================================================
//  7-node triangular element with 3-point gaussian quadrature
//
class RgTri7G3 : public RgTri7_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	RgTri7G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	Ai;
};

//=============================================================================
//  7-node triangular element with 4-point gaussian quadrature
//
class RgTri7G4 : public RgTri7_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	RgTri7G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	Ai;
};


//=============================================================================
//  7-node triangular element with 7-point gaussian quadrature
//
class RgTri7G7 : public RgTri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	RgTri7G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//  7-node triangular element with 7-point Gauss-Lobatto quadrature
//
class RgTri7GL7 : public RgTri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	RgTri7GL7();

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
class RgTri10_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 10 };

public:
	RgTri10_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_TRI10, et){}

	// initialization
	void init() override;
};

//=============================================================================
//  10-node triangular element with 7-point gaussian quadrature
//
class RgTri10G7 : public RgTri10_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	RgTri10G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};


//=============================================================================
//  10-node triangular element with 12-point gaussian quadrature
//
class RgTri10G12 : public RgTri10_
{
public:
	enum { NINT = 12 };

public:
	// constructor
	RgTri10G12();

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
class RgQuad8_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 8 };

public:
	RgQuad8_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_QUAD8, et) {}

	// initialization
	void init() override;
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 9 integration points
//
class RgQuad8G9 : public RgQuad8_
{
public:
	enum { NINT = 9 };

	// constructor
	RgQuad8G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 8 nodal integration points
//
class RgQuad8NI : public RgQuad8_
{
public:
    enum { NINT = 8 };
    
    // constructor
	RgQuad8NI();
    
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
class RgQuad9_ : public RgSurfaceElementTraits
{
public:
	enum { NELN = 9 };

public:
	RgQuad9_(int ni, ElementType et) : RgSurfaceElementTraits(ni, NELN, ET_QUAD9, et) {}

	// initialization
	void init() override;
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 integration points
//
class RgQuad9G9 : public RgQuad9_
{
public:
	enum { NINT = 9 };

	// constructor
	RgQuad9G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix	m_Ai;
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 nodal integration points
//
class RgQuad9NI : public RgQuad9_
{
public:
    enum { NINT = 9 };
    
    // constructor
	RgQuad9NI();
    
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
class RgTrussElementTraits : public RgElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	RgTrussElementTraits() : RgElementTraits(NINT, NELN, FE_ELEM_TRUSS, ET_TRUSS2, FE_TRUSS) { init(); }

	void init();
};

//=============================================================================
//          D I S C R E T E    E L E M E N T S
//
// This section defines discrete elements for 3D analysis
//=============================================================================

//=============================================================================
class RgDiscreteElementTraits : public RgElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	RgDiscreteElementTraits() : RgElementTraits(NINT, NELN, FE_ELEM_DISCRETE, ET_DISCRETE, FE_DISCRETE) { init(); }

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
class Rg2DElementTraits: public RgElementTraits
{
public:
	Rg2DElementTraits(int ni, int ne, ElementShape es, ElementType et);

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
class Rg2DTri3_ : public Rg2DElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	Rg2DTri3_(int ni, ElementType et) : Rg2DElementTraits(ni, NELN, ET_TRI3, et){}

	//! shape function at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class Rg2DTri3G1 : public Rg2DTri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	Rg2DTri3G1();

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
class Rg2DTri6_ : public Rg2DElementTraits
{
public:
	enum { NELN = 6 };

public:
	Rg2DTri6_(int ni, ElementType et) : Rg2DElementTraits(ni, NELN, ET_TRI6, et){}

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
class Rg2DTri6G3 : public Rg2DTri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	Rg2DTri6G3();

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
class Rg2DQuad4_ : public Rg2DElementTraits
{
public:
	enum { NELN = 4 };

public:
	//! constructor
	Rg2DQuad4_(int ni, ElementType et) : Rg2DElementTraits(ni, NELN, ET_QUAD4, et){}

	//! shape functions at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class Rg2DQuad4G4 : public Rg2DQuad4_
{
public:
	enum { NINT = 4 };

public:
	Rg2DQuad4G4();

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
class Rg2DQuad8_ : public Rg2DElementTraits
{
public:
	enum { NELN = 8 };

public:
	Rg2DQuad8_(int ni, ElementType et) : Rg2DElementTraits(ni, NELN, ET_QUAD8, et) {}
};