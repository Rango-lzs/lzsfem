#include "RgNLShellElement.h"



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


