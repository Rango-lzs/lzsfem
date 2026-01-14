#include "materials/SmallDeformation/RgElastoPlasticMaterialPoint.h"
#include "basicio/DumpStream.h"
#include <cassert>


namespace SmallDef {

RgElastoPlasticMaterialPoint::RgElastoPlasticMaterialPoint(RgMaterialPointData* mp) 
    : SmallDefRgMaterialPointData(mp), m_ep_bar(0.0), m_bPlastic(false), m_yield_stress(0.0)
{
    // Initialize matrices to zero
    m_e.resize(6, 1);
    m_e.zero();
    
    m_ep.resize(6, 1);
    m_ep.zero();
    
    m_ee.resize(6, 1);
    m_ee.zero();
    
    m_stress.resize(6, 1);
    m_stress.zero();
}

void RgElastoPlasticMaterialPoint::init()
{
    // Initialize base class
    SmallDefRgMaterialPointData::init();
    
    // Reset all data
    m_e.zero();
    m_ep.zero();
    m_ee.zero();
    m_stress.zero();
    m_ep_bar = 0.0;
    m_bPlastic = false;
    m_yield_stress = 0.0;
}

RgMaterialPointData* RgElastoPlasticMaterialPoint::copy()
{
    return new RgElastoPlasticMaterialPoint(*this);
}

void RgElastoPlasticMaterialPoint::serialize(DumpStream& ar)
{
    // Serialize base class
    RgMaterialPointData::serialize(ar);
    
    // Serialize plasticity data
    if (ar.isSaving()) {
        ar << m_e << m_ep << m_ee << m_stress << m_ep_bar << m_bPlastic << m_yield_stress;
    } else {
        ar >> m_e >> m_ep >> m_ee >> m_stress >> m_ep_bar >> m_bPlastic >> m_yield_stress;
    }
}

const Matrix& RgElastoPlasticMaterialPoint::getTotalStrain() const
{
    return m_e;
}

void RgElastoPlasticMaterialPoint::setTotalStrain(const Matrix& e)
{
    m_e = e;
}

const Matrix& RgElastoPlasticMaterialPoint::getPlasticStrain() const
{
    return m_ep;
}

void RgElastoPlasticMaterialPoint::setPlasticStrain(const Matrix& ep)
{
    m_ep = ep;
}

const Matrix& RgElastoPlasticMaterialPoint::getElasticStrain() const
{
    return m_ee;
}

void RgElastoPlasticMaterialPoint::calculateElasticStrain()
{
    // Elastic strain = Total strain - Plastic strain
    m_ee = m_e;
    m_ee -= m_ep;
}

double RgElastoPlasticMaterialPoint::getAccumulatedPlasticStrain() const
{
    return m_ep_bar;
}

void RgElastoPlasticMaterialPoint::setAccumulatedPlasticStrain(double ep_bar)
{
    m_ep_bar = ep_bar;
}

const Matrix& RgElastoPlasticMaterialPoint::getStress() const
{
    return m_stress;
}

void RgElastoPlasticMaterialPoint::setStress(const Matrix& stress)
{
    m_stress = stress;
}

bool RgElastoPlasticMaterialPoint::isPlastic() const
{
    return m_bPlastic;
}

void RgElastoPlasticMaterialPoint::setPlastic(bool plastic)
{
    m_bPlastic = plastic;
}

} // namespace SmallDef
