#ifndef RGHEX8GEOMNLELEMENT_H
#define RGHEX8GEOMNLELEMENT_H

#include "RgNLSolid3dElement.h"
#include <array>
#include <vector>

// Forward declarations
class RgMaterial;

/*
    RgHex8GeomNLElement
    - 8-node linear hexahedral element with geometric nonlinearity
    - Derives from RgNLSolid3dElement (nonlinear 3D solid element base)
    - Used for large deformation and large displacement analysis
    - Accounts for large deformations and rotations
    - Features:
      * Own shape function implementations (8-node linear hex)
      * 8-point Gauss quadrature (2×2×2)
      * Deformation gradient F computation
      * Green-Lagrange strain (material description)
      * Cauchy stress (spatial description)
      * Geometric (initial stress) stiffness matrix
      * Supports hyperelastic and elastoplastic materials
*/
class RgHex8GeomNLElement : public RgNLSolid3dElement
{
public:
    static constexpr int kNodeCount = 8;

    // Constructors and Destructors
    RgHex8GeomNLElement();
    explicit RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds);
    RgHex8GeomNLElement(const RgHex8GeomNLElement& other);
    RgHex8GeomNLElement& operator=(const RgHex8GeomNLElement& other);
    virtual ~RgHex8GeomNLElement();

    ElementType elementType() const override;
    ElementShape elementShape() const override;
    ElementCategory elementCategory() const override;

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

protected:
    // Compute B matrix (strain-displacement matrix) at natural coordinates
    void computeBMatrix(const Vector3d& naturalCoord, Matrix& B);
};



#endif // RGHEX8GEOMNLELEMENT_H