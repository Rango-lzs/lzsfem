#pragma once

#include "materials/RgMaterial.h"


namespace LargeDef {

/// Abstract base class for finite strain materials
class RgLargeDefMaterial : public RgMaterial {
public:
    /// Compute constitutive response for finite strain materials
    /// Uses Green-Lagrange strain (E) stored in RgMaterialPointData
    /// Constitutive matrix D is returned as an output parameter
    virtual void computeConstitutive(RgMaterialPoint* mp, Matrix& D) override = 0;
};

} // namespace LargeDef
