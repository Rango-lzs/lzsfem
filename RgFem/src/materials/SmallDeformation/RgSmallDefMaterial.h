#pragma once

#include "materials/RgMaterial.h"


namespace SmallDef {

/// Abstract base class for small strain materials
class RgSmallDefMaterial : public RgMaterial {
public:
    /// Compute constitutive response for small strain materials
    /// Uses small strain (e) stored in RgMaterialPointData
    /// Constitutive matrix D is returned as an output parameter
    virtual void computeConstitutive(RgMaterialPoint* mp, Matrix& D) override = 0;
};

} // namespace SmallDef
