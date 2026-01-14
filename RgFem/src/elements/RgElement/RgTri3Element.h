#pragma once

#include "RgLinearSolid2dElement.h"
#include <array>
#include <vector>



//! RgTri3Element - 3-node linear triangular 2D element
//! 2D solid element with linear shape functions
//! Node numbering (counterclockwise):
//!   0: (0, 0)
//!   1: (1, 0)
//!   2: (0, 1)
class FEM_EXPORT RgTri3Element : public RgLinearSolid2dElement
{
public:
    static constexpr int kNodeCount = 3;
    static constexpr int kNumEdges = 3;

    // Constructors and Destructors
    RgTri3Element();
    explicit RgTri3Element(const std::array<int, kNodeCount>& nodeIds);
    RgTri3Element(const RgTri3Element& other);
    RgTri3Element& operator=(const RgTri3Element& other);
    virtual ~RgTri3Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementCategory() const;


    // Matrix calculations
    virtual void calculateStiffnessMatrix(Matrix& K) override;
    virtual void calculateMassMatrix(Matrix& M) override;
    virtual void calculateInternalForceVector(RgVector& F) override;

};