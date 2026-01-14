#include "RgBeamElement.h"
#include "basicio/DumpStream.h"

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgBeamElement::RgBeamElement()
    : RgStructureElement()
{
}

RgBeamElement::RgBeamElement(const RgBeamElement& other)
    : RgStructureElement(other)
{
}

RgBeamElement::~RgBeamElement()
{
}

RgBeamElement& RgBeamElement::operator=(const RgBeamElement& other)
{
    if (this != &other) {
        RgStructureElement::operator=(other);
    }
    return *this;
}
