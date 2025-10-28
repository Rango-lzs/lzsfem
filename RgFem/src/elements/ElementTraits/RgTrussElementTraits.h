//=============================================================================
//          T R U S S    E L E M E N T S
//
// This section defines truss elements for 3D analysis
//=============================================================================

//=============================================================================
class FETrussElementTraits : public FEElementTraits
{
public:
    enum
    {
        NINT = 1
    };
    enum
    {
        NELN = 2
    };

public:
    FETrussElementTraits()
        : FEElementTraits(NINT, NELN, FE_ELEM_TRUSS, ET_TRUSS2, FE_TRUSS)
    {
        init();
    }

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
    enum
    {
        NINT = 1
    };
    enum
    {
        NELN = 2
    };

public:
    FEDiscreteElementTraits()
        : FEElementTraits(NINT, NELN, FE_ELEM_DISCRETE, ET_DISCRETE, FE_DISCRETE)
    {
        init();
    }

    void init()
    {
    }
};


