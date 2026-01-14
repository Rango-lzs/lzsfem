#pragma once

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>



//! RgPyramid5Element - 5-node linear pyramid element
//! 3D solid element with quadrilateral base and apex
class FEM_EXPORT RgPyramid5Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 5;

    // Constructors and Destructors
    RgPyramid5Element();
    explicit RgPyramid5Element(const std::array<int, kNodeCount>& nodeIds);
    RgPyramid5Element(const RgPyramid5Element& other);
    RgPyramid5Element& operator=(const RgPyramid5Element& other);
    virtual ~RgPyramid5Element();

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


