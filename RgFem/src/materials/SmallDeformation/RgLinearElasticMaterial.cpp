#include "materials/SmallDeformation/RgLinearElasticMaterial.h"
#include "materials/MaterialPointData.h"
#include "datastructure/mathalg.h"

namespace RgFem {
namespace SmallDef {

RgLinearElastic::RgLinearElastic(double E, double nu)
    : m_E(E), m_nu(nu)
{
    computeElasticMatrix();
}

FEMaterialPointData* RgLinearElastic::createMaterialPointData() const
{
    return new MaterialPointData();
}

void RgLinearElastic::computeConstitutive(MaterialPointData* mp, Matrix& D)
{
    // For small strain, directly return the elastic matrix
    // In this implementation, we're returning the precomputed elastic matrix
    D = m_Ce;
    
    // Also compute stress based on strain
    // sigma = D * e (in Voigt notation)
    // For simplicity, we're just setting the matrix output here
}

void RgLinearElastic::commitState(MaterialPointData* mp)
{
    mp->commit();
}

void RgLinearElastic::revertState(MaterialPointData* mp)
{
    mp->revert();
}

std::string RgLinearElastic::getName() const
{
    return "RgLinearElastic";
}

void RgLinearElastic::computeElasticMatrix()
{
    // Compute the elastic matrix for isotropic linear elasticity
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
} // namespace RgFem