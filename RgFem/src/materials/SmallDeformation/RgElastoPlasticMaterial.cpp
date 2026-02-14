#include "materials/SmallDeformation/RgElastoPlasticMaterial.h"
#include "datastructure/mathalg.h"
#include <cmath>


namespace SmallDef {

RgElastoPlastic::RgElastoPlastic(double E, double nu, double sy, double H)
    : m_E(E), m_nu(nu), m_sy(sy), m_H(H)
{
    computeElasticMatrix();
}

RgMaterialPointData* RgElastoPlastic::createRgMaterialPointData() const
{
    return new RgElastoPlasticMaterialPoint();
}

void RgElastoPlastic::computeConstitutive(RgMaterialPoint* mp, Matrix& D)
{
    // Cast to our specific material point type
    RgElastoPlasticMaterialPoint* ep_pt = dynamic_cast<RgElastoPlasticMaterialPoint*>(mp);
    if (!ep_pt) return;
    
    // Check if in plastic state to determine appropriate tangent stiffness
    if (ep_pt->isPlastic()) {
        // Return algorithmic consistent tangent for plastic state
        // For radial return algorithm with von Mises plasticity, the consistent tangent
        // is a combination of elastic modulus and plastic flow direction
        
        // Get current stress state
        const Matrix& stress = ep_pt->getStress();
        
        // Calculate deviatoric stress
        Matrix deviatoricStress = calculateDeviatoricStress(stress);
        
        // Calculate von Mises stress
        double vonMisesStress = calculateVonMisesStress(deviatoricStress);
        
        // Shear modulus
        double G = m_E / (2.0 * (1.0 + m_nu));
        
        // Calculate plastic factor
        double plasticFactor = 0.0;
        if (vonMisesStress > 0.0) {
            // Consistent tangent calculation for radial return
            double H_bar = m_H + 3.0 * G; // Combined hardening and shear modulus
            plasticFactor = (3.0 * G * (1.0 - m_sy / vonMisesStress)) / H_bar;
        }
        
        // Start with elastic modulus
        D = m_Ce;
        
        // Modify for plastic behavior
        if (vonMisesStress > 0.0) {
            // Deviatoric projection tensor
            Matrix I_dev(6, 6);
            I_dev.zero();
            
            // Identity part
            I_dev[0][0] = 1.0; I_dev[1][1] = 1.0; I_dev[2][2] = 1.0;
            I_dev[3][3] = 1.0; I_dev[4][4] = 1.0; I_dev[5][5] = 1.0;
            
            // Hydrostatic part subtraction (1/3 * δ_ij * δ_kl)
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    I_dev[i][j] -= 1.0/3.0;
                }
            }
            
            // Outer product of deviatoric stress directions
            Matrix stressNorm(6, 1);
            for (int i = 0; i < 6; i++) {
                stressNorm[i][0] = deviatoricStress[i][0] / vonMisesStress;
            }
            
            // Subtract plastic correction term
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    D[i][j] -= plasticFactor * vonMisesStress * I_dev[i][j];
                    // Additional term for stress direction dependence
                    D[i][j] -= 3.0 * G * plasticFactor * stressNorm[i][0] * stressNorm[j][0];
                }
            }
        }
    } else {
        // Elastic step - return elastic matrix
        D = m_Ce;
    }
}

void RgElastoPlastic::commitState(RgMaterialPoint* mp)
{
    
}

void RgElastoPlastic::revertState(RgMaterialPoint* mp)
{
    //mp->revert();
}

std::string RgElastoPlastic::getName() const
{
    return "RgElastoPlastic";
}

void RgElastoPlastic::updateState(RgMaterialPoint* mp, const Matrix& strainIncrement)
{
    // Cast to our specific material point type
    RgElastoPlasticMaterialPoint* ep_pt = dynamic_cast<RgElastoPlasticMaterialPoint*>(mp);
    if (!ep_pt) return;
    
    // Store the total strain
    Matrix totalStrain = ep_pt->getTotalStrain();
    totalStrain += strainIncrement;
    ep_pt->setTotalStrain(totalStrain);
    
    // Calculate trial stress
    Matrix trialStress = calculateTrialStress(ep_pt, strainIncrement);
    
    // Check for yielding
    double kappa = ep_pt->getAccumulatedPlasticStrain() * m_H; // Hardening contribution
    if (isYielding(trialStress, kappa)) {
        // Plastic step - perform return mapping
        performRadialReturn(ep_pt, trialStress);
        ep_pt->setPlastic(true);
    } else {
        // Elastic step - trial stress is the final stress
        ep_pt->setStress(trialStress);
        ep_pt->setPlastic(false);
        
        // Update elastic strain
        ep_pt->calculateElasticStrain();
    }
}

void RgElastoPlastic::computeElasticMatrix()
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

Matrix RgElastoPlastic::calculateTrialStress(RgElastoPlasticMaterialPoint* ep_pt, const Matrix& strainIncrement) const
{
    // Get current elastic strain
    Matrix elasticStrain = ep_pt->getElasticStrain();
    
    // Add strain increment to get new elastic trial strain
    elasticStrain += strainIncrement;
    
    // Calculate trial stress using Hooke's law: σ_trial = C : ε_elastic_trial
    Matrix trialStress(6, 1);
    trialStress.zero();
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            trialStress[i][0] += m_Ce[i][j] * elasticStrain[j][0];
        }
    }
    
    return trialStress;
}

bool RgElastoPlastic::isYielding(const Matrix& stress, double kappa) const
{
    // von Mises yield criterion with hardening
    double vonMisesStress = calculateVonMisesStress(stress);
    double yieldStress = m_sy + kappa; // Initial yield stress + hardening
    
    return vonMisesStress >= yieldStress;
}

void RgElastoPlastic::performRadialReturn(RgElastoPlasticMaterialPoint* ep_pt, const Matrix& trialStress) const
{
    // Calculate deviatoric trial stress
    Matrix deviatoricStress = calculateDeviatoricStress(trialStress);
    
    // Calculate von Mises trial stress
    double vonMisesStress = calculateVonMisesStress(deviatoricStress);
    
    // Current accumulated plastic strain
    double ep_bar = ep_pt->getAccumulatedPlasticStrain();
    
    // Calculate plastic multiplier using Newton-Raphson iteration
    double deltaGamma = 0.0;
    double residual = 0.0;
    const double tolerance = 1e-12;
    const int maxIterations = 20;
    
    // Shear modulus
    double G = m_E / (2.0 * (1.0 + m_nu));
    
    // Perform Newton-Raphson iterations to solve for plastic multiplier
    for (int iter = 0; iter < maxIterations; iter++) {
        double yieldFunction = vonMisesStress - 3.0 * G * deltaGamma - m_sy - m_H * ep_bar;
        residual = fabs(yieldFunction);
        
        if (residual < tolerance) break;
        
        // Derivative of yield function w.r.t. deltaGamma
        double derivative = -3.0 * G - m_H;
        
        // Update deltaGamma
        deltaGamma -= yieldFunction / derivative;
    }
    
    // Calculate final stress using radial return
    Matrix finalStress(6, 1);
    if (vonMisesStress > 0.0) {
        double scale = 1.0 - (3.0 * G * deltaGamma) / vonMisesStress;
        for (int i = 0; i < 6; i++) {
            finalStress[i][0] = scale * deviatoricStress[i][0];
        }
        // Add back hydrostatic stress
        double hydrostatic = (trialStress[0][0] + trialStress[1][0] + trialStress[2][0]) / 3.0;
        finalStress[0][0] += hydrostatic;
        finalStress[1][0] += hydrostatic;
        finalStress[2][0] += hydrostatic;
    } else {
        finalStress = trialStress;
    }
    
    // Update material point data
    ep_pt->setStress(finalStress);
    
    // Update plastic strain increment
    Matrix plasticStrainIncrement(6, 1);
    if (vonMisesStress > 0.0) {
        double scale = deltaGamma * 3.0 / vonMisesStress;
        for (int i = 0; i < 6; i++) {
            plasticStrainIncrement[i][0] = scale * deviatoricStress[i][0];
        }
    }
    
    Matrix plasticStrain = ep_pt->getPlasticStrain();
    plasticStrain += plasticStrainIncrement;
    ep_pt->setPlasticStrain(plasticStrain);
    
    // Update accumulated plastic strain
    double deltaEpBar = deltaGamma;
    ep_pt->setAccumulatedPlasticStrain(ep_bar + deltaEpBar);
    
    // Update elastic strain
    ep_pt->calculateElasticStrain();
}

double RgElastoPlastic::calculateVonMisesStress(const Matrix& stress) const
{
    // Calculate von Mises equivalent stress
    double sxx = stress[0][0];
    double syy = stress[1][0];
    double szz = stress[2][0];
    double sxy = stress[3][0];
    double syz = stress[4][0];
    double sxz = stress[5][0];
    
    double vonMises = sqrt(0.5 * (pow(sxx-syy, 2) + pow(syy-szz, 2) + pow(szz-sxx, 2)) +
                          3.0 * (sxy*sxy + syz*syz + sxz*sxz));
    
    return vonMises;
}

Matrix RgElastoPlastic::calculateDeviatoricStress(const Matrix& stress) const
{
    // Calculate hydrostatic stress
    double pressure = (stress[0][0] + stress[1][0] + stress[2][0]) / 3.0;
    
    // Calculate deviatoric stress
    Matrix deviatoric(6, 1);
    deviatoric[0][0] = stress[0][0] - pressure;
    deviatoric[1][0] = stress[1][0] - pressure;
    deviatoric[2][0] = stress[2][0] - pressure;
    deviatoric[3][0] = stress[3][0];
    deviatoric[4][0] = stress[4][0];
    deviatoric[5][0] = stress[5][0];
    
    return deviatoric;
}

} // namespace SmallDef
