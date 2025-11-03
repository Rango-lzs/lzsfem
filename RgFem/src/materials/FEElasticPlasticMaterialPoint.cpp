#include "FEElasticPlasticMaterialPoint.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEElasticPlasticMaterialPoint::FEElasticPlasticMaterialPoint(FEMaterialPointData* mp) : FEElasticMaterialPoint(mp)
{
    m_ep.zero();
    m_ee.zero();
    m_ep_mag = 0.0;
    m_bPlastic = false;
    m_yield_stress = 0.0;
}

//-----------------------------------------------------------------------------
void FEElasticPlasticMaterialPoint::Init()
{
    // Initialize base class
    FEElasticMaterialPoint::Init();
    
    // Initialize plasticity data
    m_ep.zero();
    m_ee.zero();
    m_ep_mag = 0.0;
    m_bPlastic = false;
    m_yield_stress = 0.0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEElasticPlasticMaterialPoint::Copy()
{
    FEElasticPlasticMaterialPoint* pt = new FEElasticPlasticMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEElasticPlasticMaterialPoint::Serialize(DumpStream& ar)
{
    // Serialize base class
    FEElasticMaterialPoint::Serialize(ar);
    
    // Serialize plasticity data
    ar & m_ep & m_ee & m_ep_mag & m_bPlastic & m_yield_stress;
}

//-----------------------------------------------------------------------------
Matrix3ds FEElasticPlasticMaterialPoint::PlasticStrain() const
{
    return m_ep;
}

//-----------------------------------------------------------------------------
Matrix3ds FEElasticPlasticMaterialPoint::ElasticStrain() const
{
    return m_ee;
}

//-----------------------------------------------------------------------------
void FEElasticPlasticMaterialPoint::SetPlasticStrain(const Matrix3ds& ep)
{
    m_ep = ep;
}

//-----------------------------------------------------------------------------
void FEElasticPlasticMaterialPoint::SetElasticStrain(const Matrix3ds& ee)
{
    m_ee = ee;
}

//-----------------------------------------------------------------------------
double FEElasticPlasticMaterialPoint::PlasticStrainMagnitude() const
{
    return m_ep_mag;
}

//-----------------------------------------------------------------------------
void FEElasticPlasticMaterialPoint::SetPlasticStrainMagnitude(double ep_mag)
{
    m_ep_mag = ep_mag;
}

//-----------------------------------------------------------------------------
bool FEElasticPlasticMaterialPoint::IsPlastic() const
{
    return m_bPlastic;
}