#include "RgHex20Element.h"
#include "elements/ElementShape/RgElementShape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/RgMaterialPoint.h"
#include "materials/RgMaterial.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include "femcore/RgElementTraitsStore.h"
#include "elements/NaturalCoord.h"
#include <cmath>
#include <array>


// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex20Element::RgHex20Element()
    : RgLinearSolid3dElement()
{
    initTraits();
}

RgHex20Element::RgHex20Element(const std::array<int, kNodeCount>& nodeIds)
    : RgLinearSolid3dElement()
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
    initTraits();
}

RgHex20Element::RgHex20Element(const RgHex20Element& other)
    : RgLinearSolid3dElement(other)
{
}

RgHex20Element& RgHex20Element::operator=(const RgHex20Element& other)
{
    if (this != &other) {
        RgLinearSolid3dElement::operator=(other);
    }
    return *this;
}

RgHex20Element::~RgHex20Element()
{
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgHex20Element::elementType() const
{
    return ElementType::FE_HEX20G8;
}

ElementShape RgHex20Element::elementShape() const
{
    return ElementShape::ET_HEX20;
}

ElementCategory RgHex20Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

// ============================================================================
// B-Matrix Computation
// ============================================================================

void RgHex20Element::computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B)
{
}

// ============================================================================
// Serialization
// ============================================================================

void RgHex20Element::Serialize(DumpStream& ar)
{
    RgLinearSolid3dElement::Serialize(ar);
}
