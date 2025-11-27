//=============================================================================
//                      2 D   E L E M E N T S
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
// This class defines the traits for 2D elements and serves as a
// base class for the specific 2D element formulations.
class FE2DElementTraits : public RgElementTraits
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
    enum
    {
        NELN = 3
    };

public:
    //! constructor
    FE2DTri3_(int ni, ElementType et)
        : FE2DElementTraits(ni, NELN, ET_TRI3, et)
    {
    }

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
    enum
    {
        NINT = 1
    };

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
    enum
    {
        NELN = 6
    };

public:
    FE2DTri6_(int ni, ElementType et)
        : FE2DElementTraits(ni, NELN, ET_TRI6, et)
    {
    }

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
    enum
    {
        NINT = 3
    };

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
    enum
    {
        NELN = 4
    };

public:
    //! constructor
    FE2DQuad4_(int ni, ElementType et)
        : FE2DElementTraits(ni, NELN, ET_QUAD4, et)
    {
    }

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
    enum
    {
        NINT = 4
    };

public:
    FE2DQuad4G4();

    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;

protected:
    Matrix m_Hi;  //!< inverse of H; useful for projection integr. point data to nodal data
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
    enum
    {
        NELN = 8
    };

public:
    FE2DQuad8_(int ni, ElementType et)
        : FE2DElementTraits(ni, NELN, ET_QUAD8, et)
    {
    }

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
    enum
    {
        NINT = 9
    };

    // constructor
    FE2DQuad8G9();

    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;

private:
    Matrix m_Ai;
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
    enum
    {
        NELN = 9
    };

public:
    FE2DQuad9_(int ni, ElementType et)
        : FE2DElementTraits(ni, NELN, ET_QUAD9, et)
    {
    }

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
    enum
    {
        NINT = 9
    };

    // constructor
    FE2DQuad9G9();

    // project integration point data to nodes
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
    std::vector<double> gr;  //!< integration point coordinates
    std::vector<double> gw;  //!< integration point weights

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
    enum
    {
        NELN = 2
    };

public:
    //! constructor
    FELine2_(int ni, ElementType et)
        : FELineElementTraits(ni, NELN, ET_LINE2, et)
    {
    }

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
    enum
    {
        NINT = 1
    };

public:
    //! constructor
    FELine2G1();

    //! project integration point data to nodes
    void project_to_nodes(double* ai, double* ao) const override;
};