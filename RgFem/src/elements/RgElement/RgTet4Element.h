#ifndef RGTET4ELEMENT_H
#define RGTET4ELEMENT_H

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>



// forward declarations
class Matrix3ds;

/*
    RgTet4Element
    - 4-node linear tetrahedral solid element
    - Isoparametric element with linear shape functions
    - Derives from RgLinearSolid3dElement
    - Features:
      * Constant strain element (C0 element)
      * Full integration (1 Gauss point)
      * Supports small strain analysis
      * Suitable for mesh generation and computational efficiency
*/
class RgTet4Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 4;

    // Constructors and Destructors
    RgTet4Element();
    explicit RgTet4Element(const std::array<int, kNodeCount>& nodeIds);
    RgTet4Element(const RgTet4Element& other);
    RgTet4Element& operator=(const RgTet4Element& other);
    virtual ~RgTet4Element();

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



#endif // RGTET4ELEMENT_H
