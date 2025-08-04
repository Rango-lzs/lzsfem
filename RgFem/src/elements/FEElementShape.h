#pragma once
#include "femcore/fecore_enum.h"
#include "RgElemTypeDefine.h"

//=============================================================================
// Base class for defining element shape classes.
// The element shape classes are used for evaluating shape functions and their derivatives.
class FEElementShape
{
public:
    FEElementShape(ElementShape eshape, int nodes);
    virtual ~FEElementShape();

    ElementShape shape() const
    {
        return m_shape;
    }

    int nodes() const
    {
        return m_nodes;
    }

private:
    ElementShape m_shape;
    int m_nodes;
};
