#include "RgNLShellElement.h"

namespace RgFem {

RgNLShellElement::RgNLShellElement()
    : RgElement()
{
}

RgNLShellElement::RgNLShellElement(const RgNLShellElement& other)
    : RgElement(other)
{
}

RgNLShellElement& RgNLShellElement::operator=(const RgNLShellElement& other)
{
    if (this != &other) {
        RgElement::operator=(other);
    }
    return *this;
}

RgNLShellElement::~RgNLShellElement()
{
}

} // namespace RgFem
