
#include "FEElementShape.h"
#include "elements/ElementShape/RgElementShape.h"

RgElementShape::RgElementShape(ElementShape shape, int nodes) : mShpType(shape), mNodes(nodes)
{
}

FEElementShape::~FEElementShape() 
{
}
