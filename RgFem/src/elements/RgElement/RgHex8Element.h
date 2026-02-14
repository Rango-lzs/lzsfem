#ifndef RGHEX8ELEMENT_H
#define RGHEX8ELEMENT_H

#include "datastructure/Vector3d.h"
#include "RgLinearSolid3dElement.h"

#include <array>
#include <string>
#include <vector>

/*
    RgHex8Element
    - 8-node linear hexahedral solid element (Q1 element)
    - Isoparametric element with trilinear shape functions
    - Derivatives from RgSolid3dElement
    - Features:
      * Reduced integration (1 Gauss point) or full integration (8 Gauss points)
      * Supports small strain and geometrically nonlinear analysis
      * Suitable for 3D continuum mechanics problems
*/
class RgHex8Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 8;

    // Constructors and Destructors
    RgHex8Element();
    explicit RgHex8Element(bool fullInt);
    RgHex8Element(const RgHex8Element& other);
    RgHex8Element& operator=(const RgHex8Element& other);
    virtual ~RgHex8Element();

    ElementType elementType() const override;
    ElementShape elementShape() const override;
    ElementCategory elementCategory() const override;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

protected:
    // Compute B matrix (strain-displacement matrix) at natural coordinates
    void computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B) override;
};

#endif  // RGHEX8ELEMENT_H