#pragma once
#include "datastructure/Matrix3d.h"
#include "datastructure/tens4d.h"
#include "materials/LargeDeformation/RgLargeDefMaterial.h"
#include "materials/RgMaterialPointData.h"


namespace LargeDef
{
    /// Neo-Hookean hyperelastic material model for finite strain analysis
    class RgNeoHookean : public RgLargeDefMaterial
    {
    private:
        double m_mu;      ///< Shear modulus
        double m_lambda;  ///< Lame's first parameter

    public:
        RgNeoHookean(double mu, double lambda);

        /// Create material point data
        RgMaterialPointData* createRgMaterialPointData() const override;

        /// Compute constitutive response
        void computeConstitutive(RgMaterialPoint* mp, Matrix& D) override;

        /// Commit state variables
        void commitState(RgMaterialPoint* mp) override;

        /// Revert state variables
        void revertState(RgMaterialPoint* mp) override;

        /// Get material name
        std::string getName() const override;

    private:
        /// Compute PK2 stress from deformation gradient
        Matrix3d computePK2FromF(const Matrix3d& F, double J) const;
    };

}  // namespace LargeDef
