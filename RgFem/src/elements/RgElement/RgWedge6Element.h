#pragma once

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>



//! RgWedge6Element - 6-node linear wedge (triangular prism) element
//! 3D solid element with triangular cross-section
class FEM_EXPORT RgWedge6Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 6;

    // Constructors and Destructors
    RgWedge6Element();
    explicit RgWedge6Element(const std::array<int, kNodeCount>& nodeIds);
    RgWedge6Element(const RgWedge6Element& other);
    RgWedge6Element& operator=(const RgWedge6Element& other);
    virtual ~RgWedge6Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const override;
    virtual ElementCategory elementCategory() const override;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;
   
protected:
    // Compute B matrix (strain-displacement matrix) at natural coordinates
    void computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B);
};


