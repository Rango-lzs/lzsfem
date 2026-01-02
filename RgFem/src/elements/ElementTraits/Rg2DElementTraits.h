//=============================================================================
//                      2 D   E L E M E N T S
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

#include "femcore/fem_export.h"
#include "datastructure/Matrix.h"
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix3d.h"
#include "elements/RgElemTypeDefine.h"
#include "elements/NaturalCoord.h"
#include "elements/RgGaussPoint.h"
#include "RgElementTraits.h"
#include <vector>

namespace RgFem{
    class NaturalCoord;
}

// base class for the specific 2D element formulations.
class FEM_EXPORT Rg2DElementTraits : public RgElementTraits
{
public:
    Rg2DElementTraits(int ni, int ne, ElementShape es, ElementType et);

    // initialization
    void init() override;

    // values of shape functions with size N
    virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord);

    // values of shape function derivatives with size 3,N (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord);

    // values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) {
        // Default implementation returns empty vector
        return std::vector<std::vector<double>>();
    }

    // legacy interface for shape functions
    virtual void shape(double* H, double r, double s) = 0;
    
    // legacy interface for shape function derivatives
    virtual void shape_deriv(double* Gr, double* Gs, double r, double s) = 0;

    // legacy interface for second derivatives
    virtual void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) = 0;

public:
    // gauss-points
    std::vector<RgFem::RgGaussPoint> gaussPoints;

    // local derivatives of shape functions at gauss points
    Matrix Gr, Gs;

    // local second derivatives of shape functions at gauss points
    Matrix Grr, Gsr, Grs, Gss;
};

//=============================================================================
//   Rg2DTri3
//
//=============================================================================

class FEM_EXPORT Rg2DTri3_ : public Rg2DElementTraits
{
public:
    enum
    {
        NELN = 3
    };

public:
    Rg2DTri3_(int ni, ElementType et);

    // Implementation of base class interface
    // std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    // std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    
    // Legacy interface implementation
    void shape(double* H, double r, double s) override;
    void shape_deriv(double* Gr, double* Gs, double r, double s) override;
    void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) override;
};

class FEM_EXPORT Rg2DTri3G1 : public Rg2DTri3_
{
public:
    enum
    {
        NINT = 1
    };

public:
    Rg2DTri3G1();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//   Rg2DTri6
//
//=============================================================================

class FEM_EXPORT Rg2DTri6_ : public Rg2DElementTraits
{
public:
    enum
    {
        NELN = 6
    };

public:
    Rg2DTri6_(int ni, ElementType et);

    // Implementation of base class interface
    // std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    // std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    
    // Legacy interface implementation
    void shape(double* H, double r, double s) override;
    void shape_deriv(double* Gr, double* Gs, double r, double s) override;
    void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) override;
};

class FEM_EXPORT Rg2DTri6G3 : public Rg2DTri6_
{
public:
    enum
    {
        NINT = 3
    };

public:
    Rg2DTri6G3();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//  4-node quadrilateral elements
//
class FEM_EXPORT Rg2DQuad4_ : public Rg2DElementTraits
{
public:
    enum { NELN = 4 };
    
public:
    Rg2DQuad4_(int ni, ElementType et);
    
    // Implementation of base class interface
    // std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    // std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    
    // Legacy interface implementation
    void shape(double* H, double r, double s) override;
    void shape_deriv(double* Gr, double* Gs, double r, double s) override;
    void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) override;
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature
//
class FEM_EXPORT Rg2DQuad4G4 : public Rg2DQuad4_
{
public:
    enum { NINT = 4 };

public:
    Rg2DQuad4G4();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;

private:
    Matrix m_Ai;
};

//=============================================================================
// 8-node quadrilateral elements
//
class FEM_EXPORT Rg2DQuad8_ : public Rg2DElementTraits
{
public:
    enum { NELN = 8 };

public:
    Rg2DQuad8_(int ni, ElementType et);

    // Implementation of base class interface
    // std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    // std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;  // Using base class implementation
    
    // Legacy interface implementation
    void shape(double* H, double r, double s) override;
    void shape_deriv(double* Gr, double* Gs, double r, double s) override;
    void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) override;
};

//=============================================================================
// 8-node quadrilateral elements with 9-point gaussian quadrature
//
class FEM_EXPORT Rg2DQuad8G9 : public Rg2DQuad8_
{
public:
    enum { NINT = 9 };

public:
    Rg2DQuad8G9();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;

private:
    Matrix m_Ai;
};

//=============================================================================
//                      L I N E   E L E M E N T S
//
// This section defines a set of element formulations used to describe edges.
//=============================================================================

//=============================================================================
class FEM_EXPORT RgLineElementTraits : public RgElementTraits
{
public:
    RgLineElementTraits(int ni, int ne, ElementShape es, ElementType et);

    // initialization
    void init() override;

    // shape functions at r
    virtual void shape(double* H, double r) = 0;

    // shape function derivatives at (r)
    virtual void shape_deriv(double* Gr, double r) = 0;

    // shape function second derivatives at (r)
    virtual void shape_deriv2(double* Grr, double r) = 0;

public:
    // gauss-points
    std::vector<RgFem::RgGaussPoint> gaussPoints;

    // local derivatives of shape functions at gauss points
    Matrix Gr;

    // local second derivatives of shape functions at gauss points
    Matrix Grr;
};

//=============================================================================
//
//   RgLine2_
//
//=============================================================================

//=============================================================================
//! Base class for two-point lines
class FEM_EXPORT RgLine2_ : public RgLineElementTraits
{
public:
    enum
    {
        NELN = 2
    };

public:
    RgLine2_(int ni, ElementType et);

    //! shape function at (r)
    void shape(double* H, double r) override;

    //! shape function derivatives at (r)
    void shape_deriv(double* Gr, double r) override;

    //! shape function derivatives at (r)
    void shape_deriv2(double* Grr, double r) override;
};

//=============================================================================
class FEM_EXPORT RgLine2G1 : public RgLine2_
{
public:
    enum
    {
        NINT = 1
    };

public:
    RgLine2G1();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;

private:
    Matrix m_Ai;
};