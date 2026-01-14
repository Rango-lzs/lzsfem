#include "RgQuad4Element.h"
#include <cmath>



RgQuad4Element::RgQuad4Element()
{
}

RgQuad4Element::RgQuad4Element(const std::array<int, kNodeCount>& nodeIds)
{
    
}

RgQuad4Element::RgQuad4Element(const RgQuad4Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgQuad4Element& RgQuad4Element::operator=(const RgQuad4Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgQuad4Element::~RgQuad4Element()
{
}

ElementType RgQuad4Element::elementType() const
{
    return ElementType::FE_QUAD4G4;
}

ElementShape RgQuad4Element::elementShape() const
{
    return ElementShape::ET_QUAD4;
}

ElementCategory RgQuad4Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}


void RgQuad4Element::calculateStiffnessMatrix(RgMatrix& K)
{
    // K = integral of B^T * D * B dV (using 2×2 Gauss quadrature)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 2×2 Gauss points
}

void RgQuad4Element::calculateMassMatrix(RgMatrix& M)
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgQuad4Element::calculateInternalForceVector(RgVector& F)
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    //F.zero();
    
    // Placeholder: actual implementation needed
}


