#include "RgLinearShellElement.h"

RgLinearShellElement::RgLinearShellElement()
    : RgElement()
{
}

RgLinearShellElement::RgLinearShellElement(const RgLinearShellElement& other)
    : RgElement(other)
{
}

RgLinearShellElement& RgLinearShellElement::operator=(const RgLinearShellElement& other)
{
    if (this != &other) {
        RgElement::operator=(other);
    }
    return *this;
}

RgLinearShellElement::~RgLinearShellElement()
{
}
