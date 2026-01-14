#include "RgQuad8Element.h"
#include <cmath>



RgQuad8Element::RgQuad8Element()
{
}

RgQuad8Element::RgQuad8Element(const std::array<int, kNodeCount>& nodeIds)
{
    
}

RgQuad8Element::RgQuad8Element(const RgQuad8Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgQuad8Element& RgQuad8Element::operator=(const RgQuad8Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgQuad8Element::~RgQuad8Element()
{
}

ElementType RgQuad8Element::elementType() const
{
    return ElementType::FE2D_QUAD8G9;
}

ElementShape RgQuad8Element::elementShape() const
{
    return ElementShape::ET_QUAD8;
}

ElementCategory RgQuad8Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID_2D;
}


