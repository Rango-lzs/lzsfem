#include "RgTri3Element.h"
#include <cmath>


RgTri3Element::RgTri3Element()
{
}

RgTri3Element::RgTri3Element(const std::array<int, kNodeCount>& nodeIds)
{
    
}

RgTri3Element::RgTri3Element(const RgTri3Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgTri3Element& RgTri3Element::operator=(const RgTri3Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgTri3Element::~RgTri3Element()
{
}

ElementType RgTri3Element::elementType() const
{
    return ElementType::FE_TRI3G1;
}

ElementShape RgTri3Element::elementShape() const
{
    return ElementShape::ET_TRI3;
}

ElementCategory RgTri3Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID_2D;
}

void RgTri3Element::calculateStiffnessMatrix(Matrix& K)
{
    // K = integral of B^T * D * B dV (using 1-point Gauss quadrature for linear triangle)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 1-point Gauss point
}

void RgTri3Element::calculateMassMatrix(Matrix& M)
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgTri3Element::calculateInternalForceVector(RgVector& F)
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    
    // Placeholder: actual implementation needed
}


