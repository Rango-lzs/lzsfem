#include "RgHex8GeomNLElement.h"
#include "RgNLSolid3dElement.h"
#include "RgHex8Element.h"  // For reference to linear shape functions
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>


// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8GeomNLElement::RgHex8GeomNLElement()
    : RgNLSolid3dElement()
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds)
    : RgNLSolid3dElement()
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const RgHex8GeomNLElement& other)
    : RgNLSolid3dElement(other)
{
}

RgHex8GeomNLElement::~RgHex8GeomNLElement()
{
}

ElementType RgHex8GeomNLElement::elementType() const
{
    return ElementType::FE_HEX8G8;
}

ElementShape RgHex8GeomNLElement::elementShape() const
{
    return ElementShape::ET_HEX8;
}

ElementCategory RgHex8GeomNLElement::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

void RgHex8GeomNLElement::Serialize(DumpStream& ar)
{
}

void RgHex8GeomNLElement::computeBMatrix(const Vector3d& naturalCoord, Matrix& B)
{
}

RgHex8GeomNLElement& RgHex8GeomNLElement::operator=(const RgHex8GeomNLElement& other)
{
    if (this != &other) {
        RgNLSolid3dElement::operator=(other);
    }
    return *this;
}


