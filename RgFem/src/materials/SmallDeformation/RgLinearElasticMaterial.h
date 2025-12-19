#pragma once

#include "materials/SmallDeformation/RgSmallDefMaterial.h"
#include "materials/MaterialPointData.h"
#include "datastructure/Matrix.h"
#include "datastructure/tens4d.h"

namespace RgFem {
namespace SmallDef {

/// Linear elastic material model for small strain analysis
class RgLinearElastic : public RgSmallDefMaterial {
private:
    double m_E;      ///< Young's modulus
    double m_nu;     ///< Poisson's ratio
    Matrix m_Ce;     ///< Elastic matrix in Voigt notation (6x6)
    Tensor4d m_Ce_tensor4; ///< Elastic tensor (4th order)

public:
    RgLinearElastic(double E, double nu);

    /// Create material point data
    FEMaterialPointData* createMaterialPointData() const override;

    /// Compute constitutive response
    void computeConstitutive(MaterialPointData* mp, Matrix& D) override;

    /// Commit state variables
    void commitState(MaterialPointData* mp) override;

    /// Revert state variables
    void revertState(MaterialPointData* mp) override;

    /// Get material name
    std::string getName() const override;

private:
    /// Compute elastic matrix based on material properties
    void computeElasticMatrix();
};

} // namespace SmallDef
} // namespace RgFem