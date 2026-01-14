#include "RgTri6Element.h"
#include <cmath>



RgTri6Element::RgTri6Element()
{
}

RgTri6Element::RgTri6Element(const std::array<int, kNodeCount>& nodeIds)
{
    //setNodeIds(nodeIds);
}

RgTri6Element::RgTri6Element(const RgTri6Element& other)
    : RgLinearSolid2dElement(other)
{
}

RgTri6Element& RgTri6Element::operator=(const RgTri6Element& other)
{
    if (this != &other) {
        RgLinearSolid2dElement::operator=(other);
    }
    return *this;
}

RgTri6Element::~RgTri6Element()
{
}

ElementType RgTri6Element::elementType() const
{
    return ElementType::FE_TRI6G3;
}

ElementShape RgTri6Element::elementShape() const
{
    return ElementShape::ET_TRI6;
}

ElementCategory RgTri6Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID_2D;
}


void RgTri6Element::calculateStiffnessMatrix(RgMatrix& K)
{
    // K = integral of B^T * D * B dV (using 6-point Gauss quadrature for quadratic triangle)
    int ndofs = kNodeCount * 2;  // 2D element with 2 DOF per node
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Placeholder: actual implementation needed
    // Integration loop over 6-point Gauss points
}

void RgTri6Element::calculateMassMatrix(RgMatrix& M)
{
    // M = integral of N^T * rho * N dV
    int ndofs = kNodeCount * 2;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Placeholder: actual implementation needed
}

void RgTri6Element::calculateInternalForceVector(RgVector& F)
{
    // F_int = integral of B^T * sigma dV
    int ndofs = kNodeCount * 2;
    F.resize(ndofs);
    //F.zero();
    
    // Placeholder: actual implementation needed
}


