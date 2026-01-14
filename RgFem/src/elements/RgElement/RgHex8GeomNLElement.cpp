#include "RgHex8GeomNLElement.h"
#include "RgNLSolid3dElement.h"
#include "RgHex8Element.h"  // For reference to linear shape functions
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>



// ============================================================================
// Gauss Quadrature Data (8-point for hex)
// ============================================================================

const std::array<double, 2> RgHex8GeomNLElement::gaussPoints_1D = {
    -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)
};

const std::array<double, 2> RgHex8GeomNLElement::gaussWeights_1D = {
    1.0, 1.0
};

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8GeomNLElement::RgHex8GeomNLElement()
    : RgNLSolid3dElement()
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds)
    : RgNLSolid3dElement(nodeIds)
{
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const RgHex8GeomNLElement& other)
    : RgNLSolid3dElement(other)
{
}

RgHex8GeomNLElement::~RgHex8GeomNLElement()
{
}

RgHex8GeomNLElement& RgHex8GeomNLElement::operator=(const RgHex8GeomNLElement& other)
{
    if (this != &other) {
        RgNLSolid3dElement::operator=(other);
    }
    return *this;
}


