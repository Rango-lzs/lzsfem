#pragma once
#include "femcore/fecore_enum.h"

//=============================================================================
// Base class for defining element shape classes.
// The element shape classes are used for evaluating shape functions and their derivatives.
class FEElementShape
{
public:
    FEElementShape(FE_Element_Shape eshape, int nodes);
    virtual ~FEElementShape();

    FE_Element_Shape shape() const
    {
        return m_shape;
    }
    int nodes() const
    {
        return m_nodes;
    }

private:
    FE_Element_Shape m_shape;
    int m_nodes;
};
