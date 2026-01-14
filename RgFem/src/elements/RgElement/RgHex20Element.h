#ifndef RGHEX20ELEMENT_H
#define RGHEX20ELEMENT_H

#include <array>
#include <vector>
#include "RgLinearSolid3dElement.h"



// forward declarations
class Matrix3ds;

/*
    RgHex20Element
    - 20-node serendipity hexahedral solid element (Q2 element)
    - Isoparametric element with quadratic shape functions
    - Derives from RgSolid3dElement
    - Features:
      * Biquadratic-bilinear interpolation
      * Full integration (8 Gauss points)
      * Supports small strain analysis
      * Better accuracy for curved boundaries
      * Node distribution: 8 corner + 12 mid-edge nodes
*/
class RgHex20Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 20;
    static constexpr int kNumFaces = 6;
    static constexpr int kNumEdges = 12;

    // Constructors and Destructors
    RgHex20Element();
    explicit RgHex20Element(const std::array<int, kNodeCount>& nodeIds);
    RgHex20Element(const RgHex20Element& other);
    RgHex20Element& operator=(const RgHex20Element& other);
    virtual ~RgHex20Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const override;
    virtual ElementCategory elementCategory() const override;

    virtual void Serialize(DumpStream& ar) override;

protected:
    // Compute B matrix (strain-displacement matrix)
    void computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B);
 
};


#endif // RGHEX20ELEMENT_H