#pragma once

#include "materials/RgMaterial.h"

namespace RgFem {
namespace LargeDef {

/// Abstract base class for finite strain materials
class RgLargeDefMaterial : public RgMaterial {
public:
    /// Compute constitutive response for finite strain materials
    /// Uses Green-Lagrange strain (E) stored in MaterialPointData
    /// Constitutive matrix D is returned as an output parameter
    virtual void computeConstitutive(RgMaterialPointData* mp, Matrix& D) override = 0;
};

} // namespace LargeDef
} // namespace RgFem