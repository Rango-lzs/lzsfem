#ifndef RGTET10ELEMENT_H
#define RGTET10ELEMENT_H

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>

/*
    RgTet10Element
    - 10-node quadratic tetrahedral solid element
    - Isoparametric element with quadratic shape functions
    - Derives from RgSolid3dElement
    - Features:
      * Linear strain element (C1 continuity)
      * Full integration (4 Gauss points)
      * Supports small strain analysis
      * Good accuracy for stress analysis
      * Node distribution: 4 corner + 6 mid-edge nodes
*/
class RgTet10Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 10;
    static constexpr int kNumFaces = 4;
    static constexpr int kNumEdges = 6;

    // Constructors and Destructors
    RgTet10Element();
    explicit RgTet10Element(const std::array<int, kNodeCount>& nodeIds);
    RgTet10Element(const RgTet10Element& other);
    RgTet10Element& operator=(const RgTet10Element& other);
    virtual ~RgTet10Element();

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



#endif // RGTET10ELEMENT_H
