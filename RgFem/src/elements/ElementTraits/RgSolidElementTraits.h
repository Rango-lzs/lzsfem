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
//      S O L I D   E L E M E N T 
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
//! This class defines the specific traits for solid elements and serves as
//! a base class for specific solid element formulations
//
class FEM_EXPORT FESolidElementTraits : public FEElementTraits
{
public:
	//! constructor
	FESolidElementTraits(int ni, int ne, ElementShape es, ElementType et);

	//! initialize element traits data
	void init() override;

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);

	int ShapeFunctions(int order) override;

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// element shape class
	FESolidElementShape*				m_shape;
	std::vector<FESolidElementShape*>	m_shapeP; // shape classes for different order (some orders can by null)

	// local derivatives of shape functions at gauss points
	Matrix m_Gr, m_Gs, m_Gt;
	std::vector<Matrix>	m_Gr_p;
	std::vector<Matrix>	m_Gs_p;
	std::vector<Matrix>	m_Gt_p;

	// local second derivatives of shape functions at gauss points
	Matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};


//=============================================================================
//! Base class for 8-node hexahedral elements
class FEHex8_ : public FESolidElementTraits
{
public:
	enum { NELN = 8 };

	//! initialize element traits data
	void init();

public:
	FEHex8_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_HEX8, et) {}
};

//=============================================================================
// 8-node hexahedral elements with 8-point gaussian quadrature
//
class FEHex8G8 : public FEHex8_
{
public:
	enum { NINT = 8 };

public:
	FEHex8G8();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 8-node hexahedral elements with 6-point reduced integration rule
//
class FEHex8RI : public FEHex8_
{
public:
	enum { NINT = 6 };

public:
	FEHex8RI();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 8-node hexahedral element with uniform deformation gradient

class FEHex8G1 : public FEHex8_
{
public:
	enum { NINT = 1 };

public:
	FEHex8G1();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//! Base class for 4-node linear tetrahedrons
class FETet4_ : public FESolidElementTraits
{
public:
	enum { NELN = 4 };

	//! initialize element traits data
	void init();

public:
	FETet4_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_TET4, et) {}
};

//=============================================================================
// single Gauss point integrated tet element
class FETet4G1 : public FETet4_
{
public:
	enum { NINT = 1};

public:
	FETet4G1();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 4-node tetrahedral element using a 4-node Gaussian integration rule
class FETet4G4 : public FETet4_
{
public:
	enum { NINT = 4 };

public:
	FETet4G4();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//! Base class for 5-node linear tetrahedrons
class FETet5_ : public FESolidElementTraits
{
public:
	enum { NELN = 5 };

	//! initialize element traits data
	void init();

public:
	FETet5_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_TET5, et) {}
};

//=============================================================================
// 5-node tetrahedral element using a 4-node Gaussian integration rule
class FETet5G4 : public FETet5_
{
public:
	enum { NINT = 4 };

public:
	FETet5G4();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FEPenta6
//
//=============================================================================

//=============================================================================
//! Base class for 6-node pentahedral "wedge" elements
class FEPenta6_ : public FESolidElementTraits
{
public:
	enum { NELN = 6 };

	//! initialize element traits data
	void init();

public:
	FEPenta6_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_PENTA6, et){}
};

//=============================================================================
// 6-node pentahedral elements with 6-point gaussian quadrature 
class FEPenta6G6 : public FEPenta6_
{
public:
	enum { NINT = 6 };

public:
	FEPenta6G6();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   FEPenta15
//
//=============================================================================

//=============================================================================
//! Base class for 15-node quadratic pentahedral "wedge" elements
class FEPenta15_ : public FESolidElementTraits
{
public:
    enum { NELN = 15 };
    
public:
    FEPenta15_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_PENTA15, et){}
};

//=============================================================================
// 15-node pentahedral elements with 8-point gaussian quadrature
class FEPenta15G8 : public FEPenta15_
{
public:
    enum { NINT = 8 };
    
public:
    FEPenta15G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
    Matrix m_MT;
};

//=============================================================================
// 15-node pentahedral elements with 21-point gaussian quadrature
class FEPenta15G21 : public FEPenta15_
{
public:
    enum { NINT = 21 };
    
public:
    FEPenta15G21();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
//
//   FETet10
//   
//=============================================================================

//=============================================================================
//! Base class for 10-node quadratic tetrahedral elements
class FETet10_ : public FESolidElementTraits
{
public:
	enum { NELN = 10 };

	//! initialize element traits data
	void init();

public:
	FETet10_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_TET10, et){}
};

//=============================================================================
// 10-node tetrahedral element using a 4-node Gaussian integration rule
class FETet10G1 : public FETet10_
{
public:
	enum { NINT = 1 };

public:
	FETet10G1();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix Ai;
};

//=============================================================================
// 10-node tetrahedral element using a 4-node Gaussian integration rule
class FETet10G4 : public FETet10_
{
public:
	enum { NINT = 4 };

public:
	FETet10G4();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix Ai;
};

//=============================================================================
// 10-node tetrahedral element using a 8-node Gaussian integration rule
class FETet10G8 : public FETet10_
{
public:
	enum { NINT = 8 };

public:
	FETet10G8();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
class FETet10G4RI1 : public FETet10G4, public FESRISolidElementTraits
{
public:
	FETet10G4RI1();
};

//=============================================================================
class FETet10G8RI4 : public FETet10G8, public FESRISolidElementTraits
{
public:
	FETet10G8RI4();
};

//=============================================================================
// 10-node tetrahedral element using a 11-node Gauss-Lobatto integration rule
class FETet10GL11 : public FETet10_
{
public:
	enum { NINT = 11 };

public:
	FETet10GL11();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FETet15
//   
//=============================================================================

//=============================================================================
//! Base class for 15-node quadratic tetrahedral elements
class FETet15_ : public FESolidElementTraits
{
public:
	enum { NELN = 15 };

public:
	FETet15_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_TET15, et){}
};

//=============================================================================
// 15-node tetrahedral element using a 8-node Gaussian integration rule
class FETet15G4 : public FETet15_
{
public:
	enum { NINT = 4 };

public:
	FETet15G4();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 8-node Gaussian integration rule
class FETet15G8 : public FETet15_
{
public:
	enum { NINT = 8 };

public:
	FETet15G8();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 11-point Gaussian integration rule
class FETet15G11 : public FETet15_
{
public:
	enum { NINT = 11 };

public:
	FETet15G11();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 15-point Gaussian integration rule
class FETet15G15 : public FETet15_
{
public:
	enum { NINT = 15 };

public:
	FETet15G15();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
class FETet15G15RI4 : public FETet15G15, public FESRISolidElementTraits
{
public:
	FETet15G15RI4();
};

//=============================================================================
//
//   FETet20
//   
//=============================================================================

//=============================================================================
//! Base class for 20-node cubic tetrahedral elements
class FETet20_ : public FESolidElementTraits
{
public:
	enum { NELN = 20 };

public:
	FETet20_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_TET20, et){}
};

//=============================================================================
// 20-node tetrahedral element using a 15-point Gaussian integration rule
class FETet20G15 : public FETet20_
{
public:
	enum { NINT = 15 };

public:
	FETet20G15();

	void project_to_nodes(double* ai, double* ao) const override;

private:
	Matrix N;
	Matrix Ai;
};

//=============================================================================
//
//   FEHex20
//   
//=============================================================================


//=============================================================================
//! Base class for 20-node quadratic hexahedral element
class FEHex20_ : public FESolidElementTraits
{
public:
	enum { NELN = 20 };

public:
	FEHex20_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_HEX20, et){}

	//! initialize element traits data
	void init();

	int Nodes(int order);
};

//=============================================================================
// 20-node hexahedral element using a 2x2x2 Gaussian integration rule
class FEHex20G8 : public FEHex20_
{
public:
    enum { NINT = 8 };
    
public:
    FEHex20G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    Matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
    Matrix MT;
};

//=============================================================================
// 20-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex20G27 : public FEHex20_
{
public:
	enum { NINT = 27 };

public:
	FEHex20G27();

	void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
//
//   FEHex27
//   
//=============================================================================


//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27_ : public FESolidElementTraits
{
public:
	enum { NELN = 27 };

public:
	FEHex27_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_HEX27, et){}
    
	//! initialize element traits data
	void init();

protected:
    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 27-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex27G27 : public FEHex27_
{
public:
	enum { NINT = 27 };

public:
	FEHex27G27();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 
// FEPyra5
//
//=============================================================================

//=============================================================================
//! Base class for 5-node pyramid element
class FEPyra5_ : public FESolidElementTraits
{
public:
	enum { NELN = 5 };

public:
	FEPyra5_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_PYRA5, et){}
};

//=============================================================================
// 5-node pyramid element using a 2x2x2 Gaussian integration rule
class FEPyra5G8: public FEPyra5_
{
public:
	enum { NINT = 8 };

public:
	FEPyra5G8();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	Matrix	m_Ai;
};

//=============================================================================
//
// FEPyra13
//
//=============================================================================

//=============================================================================
//! Base class for 13-node pyramid element
class FEPyra13_ : public FESolidElementTraits
{
public:
    enum { NELN = 13 };
    
public:
    FEPyra13_(int ni, ElementType et) : FESolidElementTraits(ni, NELN, ET_PYRA13, et){}
};

//=============================================================================
// 13-node pyramid element using a 2x2x2 Gaussian integration rule
class FEPyra13G8: public FEPyra13_
{
public:
    enum { NINT = 8 };
    
public:
    FEPyra13G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    Matrix    m_Ai;
};
