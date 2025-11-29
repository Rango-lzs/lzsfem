#include "RgLinearBeamElement.h"

namespace RgFem {

RgLinearBeamElement::RgLinearBeamElement()
    : RgBeamElement()
{
}

RgLinearBeamElement::RgLinearBeamElement(const RgLinearBeamElement& other)
    : RgBeamElement(other)
{
}

RgLinearBeamElement& RgLinearBeamElement::operator=(const RgLinearBeamElement& other)
{
    if (this != &other) {
        RgBeamElement::operator=(other);
    }
    return *this;
}

RgLinearBeamElement::~RgLinearBeamElement()
{
}

} // namespace RgFem
