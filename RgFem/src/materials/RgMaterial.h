#pragma once

#include "FEMaterial.h"
#include "RgMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/tens4d.h"
#include "datastructure/Matrix.h"
#include "femcore/FETimeInfo.h"
#include <vector>
#include <memory>

// Forward declarations
class RgMaterialPointData;

namespace RgFem {

// ============================================================================
// ENUMERATIONS AND TYPE DEFINITIONS
// ============================================================================

/// Strain measures supported by the material framework
enum class StrainMeasure {
    GreenLagrange,   ///< Green-Lagrange strain (E = 0.5(C-I))
    Almansi,         ///< Almansi strain (e = 0.5(I-c^-1))
    Logarithmic,     ///< Logarithmic strain (log(U))
    SmallStrain      ///< Infinitesimal strain (eps = 0.5(Grad[u] + Grad[u]^T))
};

/// Stress measures supported by the material framework
enum class StressMeasure {
    PK2,             ///< Second Piola-Kirchhoff stress
    Cauchy,          ///< Cauchy (true) stress
    Kirchhoff,       ///< Kirchhoff stress (tau = J*sigma)
    PK1              ///< First Piola-Kirchhoff stress
};

/// Update modes for finite element analysis
enum class UpdateMode {
    TotalLagrangian,    ///< Total Lagrangian formulation
    UpdatedLagrangian   ///< Updated Lagrangian formulation
};

/// Abstract base class for all material models
class RgMaterial 
{
public:
    virtual ~RgMaterial() = default;

    /// Create material point data for this material type
    virtual RgMaterialPointData* createMaterialPointData() const = 0;

    /// Compute constitutive response for given material point
    /// Material uses kinematics stored in MaterialPointData to compute stress
    /// Constitutive matrix D is returned as an output parameter
    virtual void computeConstitutive(RgMaterialPointData* mp, Matrix& D) = 0;

    /// Commit state variables at material point
    virtual void commitState(RgMaterialPointData* mp) = 0;

    /// Revert state variables at material point to last committed state
    virtual void revertState(RgMaterialPointData* mp) = 0;

    /// Get human-readable material name
    virtual std::string getName() const = 0;
};

} // namespace RgFem