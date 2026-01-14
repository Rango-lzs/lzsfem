#include "RgNLBeamElement.h"



RgNLBeamElement::RgNLBeamElement()
    : RgBeamElement()
{
}

RgNLBeamElement::RgNLBeamElement(const RgNLBeamElement& other)
    : RgBeamElement(other)
{
}

RgNLBeamElement& RgNLBeamElement::operator=(const RgNLBeamElement& other)
{
    if (this != &other) {
        RgBeamElement::operator=(other);
    }
    return *this;
}

RgNLBeamElement::~RgNLBeamElement()
{
}


