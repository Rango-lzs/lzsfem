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
class RgElement;
class RgSolidElementShape;
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

    // values of shape functions with size N
    virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord);

    // values of shape function derivatives with size 3,N (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord);

    // values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord);

    virtual void project_to_nodes(double* ai, double* ao) const;

    const RgFem::RgGaussPoint gaussPoint(int n);


public:
    // gauss-points
    std::vector<RgFem::RgGaussPoint> gaussPoints;

    // element shape class
    RgSolidElementShape* m_shape;
    // local derivatives of shape functions at gauss points
    Matrix m_Gr, m_Gs, m_Gt;
    // local second derivatives of shape functions at gauss points
    Matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};

//=============================================================================
//! Base class for 8-node hexahedral elements
class FEM_EXPORT RgHex8_ : public RgSolidElementTraits
{
public:
    enum { NELN = 8 };

    //! initialize element traits data
    void init() override;

public:
    RgHex8_(int ni, ElementType et);
};

//=============================================================================
// 8-node hexahedral elements with 8-point gaussian quadrature
//
class FEM_EXPORT RgHex8G8 : public RgHex8_
{
public:
    enum { NINT = 8 };

public:
    RgHex8G8();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi; //!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 8-node hexahedral elements with 6-point reduced integration rule
//
class FEM_EXPORT RgHex8RI : public RgHex8_
{
public:
    enum { NINT = 6 };

public:
    RgHex8RI();

    void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 8-node hexahedral element with uniform deformation gradient

class FEM_EXPORT RgHex8G1 : public RgHex8_
{
public:
    enum { NINT = 1 };

public:
    RgHex8G1();

    void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//! Base class for 4-node linear tetrahedrons
class FEM_EXPORT RgTet4_ : public RgSolidElementTraits
{
public:
    enum { NELN = 4 };

    //! initialize element traits data
    void init() override;

public:
    RgTet4_(int ni, ElementType et);
};

//=============================================================================
// single Gauss point integrated tet element
class FEM_EXPORT RgTet4G1 : public RgTet4_
{
public:
    enum { NINT = 1};

public:
    RgTet4G1();

    void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 4-node tetrahedral element using a 4-node Gaussian integration rule
class FEM_EXPORT RgTet4G4 : public RgTet4_
{
public:
    enum { NINT = 4 };

public:
    RgTet4G4();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi; //!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   RgPenta6
//
//=============================================================================

//=============================================================================
//! Base class for 6-node pentahedral "wedge" elements
class FEM_EXPORT RgPenta6_ : public RgSolidElementTraits
{
public:
    enum { NELN = 6 };

    //! initialize element traits data
    void init() override;

public:
    RgPenta6_(int ni, ElementType et);
};

//=============================================================================
// 6-node pentahedral elements with 6-point gaussian quadrature 
class FEM_EXPORT RgPenta6G6 : public RgPenta6_
{
public:
    enum { NINT = 6 };

public:
    RgPenta6G6();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi; //!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   RgPenta15
//
//=============================================================================

//=============================================================================
//! Base class for 15-node quadratic pentahedral "wedge" elements
class FEM_EXPORT RgPenta15_ : public RgSolidElementTraits
{
public:
    enum { NELN = 15 };

    //! initialize element traits data
    void init() override;

public:
    RgPenta15_(int ni, ElementType et);
};

//=============================================================================
// 15-node pentahedral elements with 21-point gaussian quadrature 
class FEM_EXPORT RgPenta15G21 : public RgPenta15_
{
public:
    enum { NINT = 21 };

public:
    RgPenta15G21();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi; //!< inverse of H; useful for projection integr. point data to nodal data 
};


//=============================================================================
//
//   RgHex20
//   
//=============================================================================


//=============================================================================
//! Base class for 20-node quadratic hexahedral element
class FEM_EXPORT RgHex20_ : public RgSolidElementTraits
{
public:
    enum { NELN = 20 };

public:
    RgHex20_(int ni, ElementType et);

    //! initialize element traits data
    void init() override;

    int Nodes(int order);
};

//=============================================================================
// 20-node hexahedral element using a 2x2x2 Gaussian integration rule
class FEM_EXPORT RgHex20G8 : public RgHex20_
{
public:
    enum { NINT = 8 };
    
public:
    RgHex20G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    Matrix Hi; //!< inverse of H; useful for projection integr. point data to nodal data
    Matrix MT;
};

//=============================================================================
// 20-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEM_EXPORT RgHex20G27 : public RgHex20_
{
public:
    enum { NINT = 27 };
    
public:
    RgHex20G27();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    Matrix Hi; //!< inverse of H; useful for projection integr. point data to nodal data
    Matrix MT;
};

//=============================================================================
//
//   RgHex27
//   
//=============================================================================

//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEM_EXPORT RgHex27_ : public RgSolidElementTraits
{
public:
    enum { NELN = 27 };

public:
    RgHex27_(int ni, ElementType et);
};

//=============================================================================
// 27-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEM_EXPORT RgHex27G27 : public RgHex27_
{
public:
    enum { NINT = 27 };

public:
    RgHex27G27();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi; //!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   RgPyra5
//
//=============================================================================

//=============================================================================
//! Base class for 5-node pyramid element
class FEM_EXPORT RgPyra5_ : public RgSolidElementTraits
{
public:
    enum { NELN = 5 };

public:
    RgPyra5_(int ni, ElementType et);
};

//=============================================================================
// 5-node pyramid element using a 2x2x2 Gaussian integration rule
class FEM_EXPORT RgPyra5G8: public RgPyra5_
{
public:
    enum { NINT = 8 };

public:
    RgPyra5G8();

    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Ai;
};

//=============================================================================
//
// RgPyra13
//
//=============================================================================

//=============================================================================
//! Base class for 13-node pyramid element
class FEM_EXPORT RgPyra13_ : public RgSolidElementTraits
{
public:
    enum { NELN = 13 };
    
public:
    RgPyra13_(int ni, ElementType et);
};

//=============================================================================
// 13-node pyramid element using a 2x2x2 Gaussian integration rule
class FEM_EXPORT RgPyra13G8: public RgPyra13_
{
public:
    enum { NINT = 8 };
    
public:
    RgPyra13G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    Matrix m_Ai;
};