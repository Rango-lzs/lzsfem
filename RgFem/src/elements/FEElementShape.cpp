
#include "FEElementShape.h"

FEElementShape::FEElementShape(ElementShape eshape, int nodes) : m_shape(eshape), m_nodes(nodes)
{
}

FEElementShape::~FEElementShape() 
{
}
