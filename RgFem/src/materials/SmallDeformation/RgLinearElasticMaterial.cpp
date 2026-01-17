#include "materials/SmallDeformation/RgLinearElasticMaterial.h"
#include "materials/SmallDeformation/RgLinearElasticMatData.h"
#include "materials/RgMaterialPointData.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/mathalg.h"


namespace SmallDef {

RgLinearElastic::RgLinearElastic(double E, double nu)
    : m_E(E), m_nu(nu)
{
    updateLameParameters();
}

void RgLinearElastic::setYoungsModulus(double E) {
    m_E = E;
    updateLameParameters();
}

void RgLinearElastic::setPoissonsRatio(double nu) {
    m_nu = nu;
    updateLameParameters();
}

double RgLinearElastic::getYoungsModulus() const {
    return m_E;
}

double RgLinearElastic::getPoissonsRatio() const {
    return m_nu;
}

void RgLinearElastic::updateLameParameters()
{
    m_lambda = (m_E * m_nu) / ((1 + m_nu) * (1 - 2 * m_nu));
    m_mu = m_E / (2 * (1 + m_nu));
}

RgMaterialPointData* RgLinearElastic::createRgMaterialPointData() const
{
    // For linear elastic materials, return an instance of the specific material point data
    return new RgLinearElasticMatData();
}

void RgLinearElastic::computeConstitutive(RgMaterialPointData* mp, Matrix& D)
{
    if (!mp) return;
    
    // For linear elasticity, we use the standard stiffness matrix
    // Calculate the 6x6 stiffness matrix for isotropic material
    D.resize(6, 6);
    D.zero();
    
    double c1 = m_E * (1 - m_nu) / ((1 + m_nu) * (1 - 2 * m_nu));  // C11
    double c12 = m_E * m_nu / ((1 + m_nu) * (1 - 2 * m_nu));       // C12
    double c44 = m_E / (2 * (1 + m_nu));                            // C44
    
    // Fill the stiffness matrix for isotropic material
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            D[i][j] = (i == j) ? c1 : c12;
        }
        D[i+3][i+3] = c44;
    }
    
    // Calculate stress from strain using the constitutive matrix
    // This would typically be done in the element implementation
    // For demonstration, we'll compute stress based on strain in mp
    if (mp) {
        // Use the strain from the material point
        // Extract the SmallDefRgMaterialPointData from the RgMaterialPointData chain
        SmallDefRgMaterialPointData* smpt = mp->ExtractData<SmallDefRgMaterialPointData>();
        if (smpt) {
            Matrix3ds strain = smpt->strain; // Using the strain field from SmallDef::SmallDefRgMaterialPointData
            Matrix3ds stress;
            
            // Compute stress = D * strain
            // This is a simplified calculation for demonstration
            stress.xx() = D[0][0] * strain.xx() + D[0][1] * strain.yy() + D[0][2] * strain.zz();
            stress.yy() = D[1][0] * strain.xx() + D[1][1] * strain.yy() + D[1][2] * strain.zz();
            stress.zz() = D[2][0] * strain.xx() + D[2][1] * strain.yy() + D[2][2] * strain.zz();
            stress.xy() = D[3][3] * strain.xy();
            stress.yz() = D[4][4] * strain.yz();
            stress.xz() = D[5][5] * strain.xz();
            
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
    return "RgLinearElastic";
}

void RgLinearElastic::computeElasticMatrix()
{
    // Compute the elastic matrix for isotropic linear elasticity in Voigt notation
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
    m_Ce[3][3] = c * (1.0 - 2.0 * m_nu) / 2.0;
    m_Ce[4][4] = c * (1.0 - 2.0 * m_nu) / 2.0;
    m_Ce[5][5] = c * (1.0 - 2.0 * m_nu) / 2.0;
}

} // namespace SmallDef
