#pragma once

#include "materials/RgMaterial.h"

namespace RgFem {
namespace SmallDef {

/// Abstract base class for small strain materials
class RgSmallDefMaterial : public RgMaterial {
public:
    /// Compute constitutive response for small strain materials
    /// Uses small strain (e) stored in MaterialPointData
    /// Constitutive matrix D is returned as an output parameter
    virtual void computeConstitutive(RgMaterialPointData* mp, Matrix& D) override = 0;
};

} // namespace SmallDef
} // namespace RgFem