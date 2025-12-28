#pragma once

#include "materials/SmallDeformation/RgSmallDefMaterial.h"
#include "materials/RgMaterialPointData.h"
#include "datastructure/Matrix.h"
#include "datastructure/tens4d.h"

namespace RgFem {
namespace SmallDef {

/// Material mode enum to distinguish between 3D, 2D plane stress and 2D plane strain
enum class MaterialMode {
    THREE_D,        ///< 3D analysis
    PLANE_STRESS,   ///< 2D plane stress analysis
    PLANE_STRAIN    ///< 2D plane strain analysis
};

/// Linear elastic material model for small strain analysis
class RgLinearElastic : public RgSmallDefMaterial {
private:
    double m_E;      ///< Young's modulus
    double m_nu;     ///< Poisson's ratio
    double m_lambda; ///< Lame's first parameter
    double m_mu;     ///< Lame's second parameter (shear modulus)
    Matrix m_Ce;     ///< Elastic matrix in Voigt notation (6x6)
    MaterialMode m_materialMode; ///< Type of material mode (3D, plane stress, or plane strain)

public:
    RgLinearElastic(double E, double nu, MaterialMode mode = MaterialMode::THREE_D);

    /// Create material point data
    RgMaterialPointData* createMaterialPointData() const override;

    /// Compute constitutive response
    void computeConstitutive(RgMaterialPointData* mp, Matrix& D) override;

    /// Commit state variables
    void commitState(RgMaterialPointData* mp) override;

    /// Revert state variables
    void revertState(RgMaterialPointData* mp) override;

    /// Get material name
    std::string getName() const override;

    /// Set Young's modulus
    void setYoungsModulus(double E);

    /// Set Poisson's ratio
    void setPoissonsRatio(double nu);

    /// Set material mode
    void setMaterialMode(MaterialMode mode);

    /// Get Young's modulus
    double getYoungsModulus() const;

    /// Get Poisson's ratio
    double getPoissonsRatio() const;

    /// Get material mode
    MaterialMode getMaterialMode() const;

private:
    /// Update Lame parameters based on E and nu
    void updateLameParameters();

    /// Compute elastic matrix based on material properties and material mode
    void computeElasticMatrix();
};

} // namespace SmallDef
} // namespace RgFem