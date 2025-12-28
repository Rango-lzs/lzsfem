#include "materials/SmallDeformation/RgLinearElasticMaterial.h"
#include "materials/SmallDeformation/RgLinearElasticMatData.h"
#include "materials/RgMaterialPointData.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/mathalg.h"

namespace RgFem {
namespace SmallDef {

RgLinearElastic::RgLinearElastic(double E, double nu, MaterialMode mode)
    : m_E(E), m_nu(nu), m_materialMode(mode)
{
    updateLameParameters();
    computeElasticMatrix();
}

void RgLinearElastic::setYoungsModulus(double E) {
    m_E = E;
    updateLameParameters();
    computeElasticMatrix();
}

void RgLinearElastic::setPoissonsRatio(double nu) {
    m_nu = nu;
    updateLameParameters();
    computeElasticMatrix();
}

void RgLinearElastic::setMaterialMode(MaterialMode mode) {
    m_materialMode = mode;
    computeElasticMatrix(); // Recompute matrix based on new material mode
}

double RgLinearElastic::getYoungsModulus() const {
    return m_E;
}

double RgLinearElastic::getPoissonsRatio() const {
    return m_nu;
}

MaterialMode RgLinearElastic::getMaterialMode() const {
    return m_materialMode;
}

void RgLinearElastic::updateLameParameters()
{
    m_lambda = (m_E * m_nu) / ((1 + m_nu) * (1 - 2 * m_nu));
    m_mu = m_E / (2 * (1 + m_nu));
}

RgMaterialPointData* RgLinearElastic::createMaterialPointData() const
{
    // For linear elastic materials, return an instance of the specific material point data
    return new RgLinearElasticMatData();
}

void RgLinearElastic::computeConstitutive(RgMaterialPointData* mp, Matrix& D)
{
    if (!mp) return;
    
    // Use the precomputed elastic matrix based on material mode
    D = m_Ce;
    
    // Calculate stress from strain using the constitutive matrix
    // This would typically be done in the element implementation
    // For demonstration, we'll compute stress based on strain in mp
    if (mp) {
        // Use the strain from the material point
        // Extract the SmallDefMaterialPointData from the RgMaterialPointData chain
        SmallDefMaterialPointData* smpt = mp->ExtractData<SmallDefMaterialPointData>();
        if (smpt) {
            Matrix3ds strain = smpt->strain; // Using the strain field from SmallDef::SmallDefMaterialPointData
            Matrix3ds stress;
            
            // Compute stress = D * strain based on material mode
            switch (m_materialMode) {
                case MaterialMode::THREE_D:
                    // 3D case: use full 6x6 matrix
                    stress.xx() = D[0][0] * strain.xx() + D[0][1] * strain.yy() + D[0][2] * strain.zz();
                    stress.yy() = D[1][0] * strain.xx() + D[1][1] * strain.yy() + D[1][2] * strain.zz();
                    stress.zz() = D[2][0] * strain.xx() + D[2][1] * strain.yy() + D[2][2] * strain.zz();
                    stress.xy() = D[3][3] * strain.xy();
                    stress.yz() = D[4][4] * strain.yz();
                    stress.xz() = D[5][5] * strain.xz();
                    break;
                    
                case MaterialMode::PLANE_STRESS:
                    // Plane stress: only in-plane stresses (xx, yy, xy)
                    // zz stress is zero, but there's zz strain
                    stress.xx() = D[0][0] * strain.xx() + D[0][1] * strain.yy();
                    stress.yy() = D[1][0] * strain.xx() + D[1][1] * strain.yy();
                    stress.xy() = D[2][2] * strain.xy();
                    stress.zz() = 0.0;  // Plane stress condition
                    stress.yz() = 0.0;
                    stress.xz() = 0.0;
                    break;
                    
                case MaterialMode::PLANE_STRAIN:
                    // Plane strain: only in-plane strains (xx, yy, xy)
                    // zz strain is zero, but there's zz stress
                    stress.xx() = D[0][0] * strain.xx() + D[0][1] * strain.yy();
                    stress.yy() = D[1][0] * strain.xx() + D[1][1] * strain.yy();
                    stress.zz() = D[2][0] * strain.xx() + D[2][1] * strain.yy();  // Non-zero due to Poisson effect
                    stress.xy() = D[3][3] * strain.xy();
                    stress.yz() = 0.0;  // Zero due to plane strain
                    stress.xz() = 0.0;  // Zero due to plane strain
                    break;
            }
            
            smpt->stress = stress; // Store result as Cauchy stress
        }
    }
}

void RgLinearElastic::commitState(RgMaterialPointData* mp)
{
    if (mp) {
        // For linear elastic, no state variables to commit
        // Just pass through to base implementation if needed
        //mp->commit();
    }
}

void RgLinearElastic::revertState(RgMaterialPointData* mp)
{
    if (mp) {
        // For linear elastic, no state variables to revert
        // Just pass through to base implementation if needed
        //mp->revert();
    }
}

std::string RgLinearElastic::getName() const
{
    std::string name = "RgLinearElastic";
    switch (m_materialMode) {
        case MaterialMode::THREE_D:
            name += " (3D)";
            break;
        case MaterialMode::PLANE_STRESS:
            name += " (Plane Stress)";
            break;
        case MaterialMode::PLANE_STRAIN:
            name += " (Plane Strain)";
            break;
    }
    return name;
}

void RgLinearElastic::computeElasticMatrix()
{
    switch (m_materialMode) {
        case MaterialMode::THREE_D:
        {
            // Compute the elastic matrix for 3D isotropic linear elasticity in Voigt notation
            const double c = m_E / ((1.0 + m_nu) * (1.0 - 2.0 * m_nu));
            
            m_Ce.resize(6, 6);
            m_Ce.zero();
            
            // Normal stresses
            m_Ce[0][0] = c * (1.0 - m_nu);
            m_Ce[1][1] = c * (1.0 - m_nu);
            m_Ce[2][2] = c * (1.0 - m_nu);
            
            m_Ce[0][1] = c * m_nu;
            m_Ce[0][2] = c * m_nu;
            m_Ce[1][0] = c * m_nu;
            m_Ce[1][2] = c * m_nu;
            m_Ce[2][0] = c * m_nu;
            m_Ce[2][1] = c * m_nu;
            
            // Shear stresses
            m_Ce[3][3] = m_E / (2.0 * (1.0 + m_nu));  // = c * (1.0 - 2.0 * m_nu) / 2.0 simplified
            m_Ce[4][4] = m_E / (2.0 * (1.0 + m_nu));
            m_Ce[5][5] = m_E / (2.0 * (1.0 + m_nu));
            break;
        }
        
        case MaterialMode::PLANE_STRESS:
        {
            // Compute the elastic matrix for 2D plane stress in Voigt notation
            // Matrix is 3x3: [sig_xx, sig_yy, sig_xy] = [D]{eps_xx, eps_yy, 2*eps_xy}
            const double factor = m_E / (1.0 - m_nu * m_nu);
            
            m_Ce.resize(3, 3);
            m_Ce.zero();
            
            // Normal stresses
            m_Ce[0][0] = factor;              // D11
            m_Ce[0][1] = factor * m_nu;       // D12
            m_Ce[1][0] = factor * m_nu;       // D21
            m_Ce[1][1] = factor;              // D22
            
            // Shear stress
            m_Ce[2][2] = m_E / (2.0 * (1.0 + m_nu));  // D33
            break;
        }
        
        case MaterialMode::PLANE_STRAIN:
        {
            // Compute the elastic matrix for 2D plane strain in Voigt notation
            // Matrix is 3x3: [sig_xx, sig_yy, sig_xy] = [D]{eps_xx, eps_yy, 2*eps_xy}
            const double c = m_E * (1.0 - m_nu) / ((1.0 + m_nu) * (1.0 - 2.0 * m_nu));
            
            m_Ce.resize(3, 3);
            m_Ce.zero();
            
            // Normal stresses
            m_Ce[0][0] = c;                             // D11
            m_Ce[0][1] = c * m_nu / (1.0 - m_nu);     // D12
            m_Ce[1][0] = c * m_nu / (1.0 - m_nu);     // D21
            m_Ce[1][1] = c;                             // D22
            
            // Shear stress
            m_Ce[2][2] = m_E / (2.0 * (1.0 + m_nu));  // D33
            break;
        }
    }
}

} // namespace SmallDef
} // namespace RgFem