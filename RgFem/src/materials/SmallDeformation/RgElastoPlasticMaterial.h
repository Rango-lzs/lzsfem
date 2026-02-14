#pragma once
#include "materials/SmallDeformation/RgSmallDefMaterial.h"
#include "materials/SmallDeformation/RgElastoPlasticMaterialPoint.h"
#include "datastructure/Matrix.h"


namespace SmallDef {

/// Elasto-plastic material model for small strain analysis
class RgElastoPlastic : public RgSmallDefMaterial {
private:
    double m_E;      ///< Young's modulus
    double m_nu;     ///< Poisson's ratio
    double m_sy;     ///< Yield stress
    double m_H;      ///< Hardening modulus
    Matrix m_Ce;     ///< Elastic matrix in Voigt notation (6x6)

public:
    RgElastoPlastic(double E, double nu, double sy, double H = 0.0);

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

    /// Update material state with strain increment
    void updateState(RgMaterialPoint* mp, const Matrix& strainIncrement);

private:
    /// Compute elastic matrix based on material properties
    void computeElasticMatrix();
    
    /// Calculate trial stress
    Matrix calculateTrialStress(RgElastoPlasticMaterialPoint* ep_pt, const Matrix& strainIncrement) const;
    
    /// Check if yielding occurs
    bool isYielding(const Matrix& stress, double kappa) const;
    
    /// Perform radial return mapping algorithm
    void performRadialReturn(RgElastoPlasticMaterialPoint* ep_pt, const Matrix& trialStress) const;
    
    /// Calculate von Mises equivalent stress
    double calculateVonMisesStress(const Matrix& stress) const;
    
    /// Calculate deviatoric stress
    Matrix calculateDeviatoricStress(const Matrix& stress) const;
};

} // namespace SmallDef
