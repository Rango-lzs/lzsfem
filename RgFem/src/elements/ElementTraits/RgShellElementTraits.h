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
#include "RgElementTraits.h"


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
class FEM_EXPORT RgShellElementTraits : public RgElementTraits
{
public:
    RgShellElementTraits(int ni, int ne, ElementShape es, ElementType et);

    void init() override;
    
    // values of shape functions with size N
    virtual std::vector<double> evalH(const NaturalCoord& coord) override;

    // values of shape function derivatives with size 3,N (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord) override;

    // values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) override ;
    
 protected:

	// local derivatives of shape functions at gauss points
	Matrix Hr, Hs;
    //Matrix Hrr, Hsr, Grs, Hss;
};

//=============================================================================
// 4-node quadrilateral elements
//
class FEM_EXPORT RgShellQuad4_ : public RgShellElementTraits
{
public:
    enum { NELN = 4 };
    
public:
    RgShellQuad4_(int ni, ElementType et);   
};

//=============================================================================
// 4-node quadrilateral elements with 4*2-point gaussian quadrature
//
class FEM_EXPORT RgShellQuad4G8 : public RgShellQuad4_
{
public:
    enum { NINT = 8 };
    
public:
    RgShellQuad4G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 4-node quadrilateral elements with 4*3-point gaussian quadrature
//
class FEM_EXPORT RgShellQuad4G12 : public RgShellQuad4_
{
public:
    enum { NINT = 12 };
    
public:
    RgShellQuad4G12();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements
//
class FEM_EXPORT RgShellTri3_ : public RgShellElementTraits
{
public:
    enum { NELN = 3 };
    
public:
    RgShellTri3_(int ni, ElementType et);
};

//=============================================================================
// 3-node triangular elements with 3*2-point gaussian quadrature
//
class FEM_EXPORT RgShellTri3G6 : public RgShellTri3_
{
public:
    enum { NINT = 6 };
    
public:
    RgShellTri3G6();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements with 3*3-point gaussian quadrature
//
class FEM_EXPORT RgShellTri3G9 : public RgShellTri3_
{
public:
    enum { NINT = 9 };
    
public:
    RgShellTri3G9();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements
//
class FEM_EXPORT RgShellQuad8_ : public RgShellElementTraits
{
public:
    enum { NELN = 8 };
    
public:
    RgShellQuad8_(int ni, ElementType et);  
};

//=============================================================================
// 8-node quadrilateral elements with 9*2-point gaussian quadrature
//
class FEM_EXPORT RgShellQuad8G18 : public RgShellQuad8_
{
public:
    enum { NINT = 18 };
    
public:
    RgShellQuad8G18();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements with 9*3-point gaussian quadrature
//
class FEM_EXPORT RgShellQuad8G27 : public RgShellQuad8_
{
public:
    enum { NINT = 27 };
    
public:
    RgShellQuad8G27();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements
//
class FEM_EXPORT RgShellTri6_ : public RgShellElementTraits
{
public:
    enum { NELN = 6 };
    
public:
    RgShellTri6_(int ni, ElementType et);       
};

//=============================================================================
// 6-node triangular elements with 7*2-point gaussian quadrature
//
class FEM_EXPORT RgShellTri6G14 : public RgShellTri6_
{
public:
    enum { NINT = 14 };
    
public:
    RgShellTri6G14();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements with 7*3-point gaussian quadrature
//
class FEM_EXPORT RgShellTri6G21 : public RgShellTri6_
{
public:
    enum { NINT = 21 };
    
public:
    RgShellTri6G21();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    Matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};