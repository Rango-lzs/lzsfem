#include "FEElasticPlasticMaterial.h"
#include "femcore/FEModel.h"
#include "datastructure/mathalg.h"

DEFINE_META_CLASS(FEElasticPlasticMaterial, FEElasticMaterial, "");

BEGIN_PARAM_DEFINE(FEElasticPlasticMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_E, FEParamDouble, "E");
	ADD_PARAMETER(m_v, FEParamDouble, "v");
	ADD_PARAMETER(m_sy, FEParamDouble, "sy");
END_PARAM_DEFINE();

FEElasticPlasticMaterial::FEElasticPlasticMaterial() : FEElasticMaterial()
{ 
	m_E = 200000;   // Default Young's modulus (Steel)
	m_v = 0.3;      // Default Poisson's ratio
	m_sy = 250;     // Default yield stress (Steel)
}

//-----------------------------------------------------------------------------
FEElasticPlasticMaterial::~FEElasticPlasticMaterial()
{ 
	
}

//-----------------------------------------------------------------------------
//! create material point data for this material
FEMaterialPointData* FEElasticPlasticMaterial::CreateMaterialPointData()
{ 
	return new FEElasticPlasticMaterialPoint;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEElasticPlasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
	FEElasticPlasticMaterialPoint& et = *pt.ExtractData<FEElasticPlasticMaterialPoint>();
	
	// Get elastic strain
	Matrix3ds e = et.ElasticStrain();
	
	// Calculate elastic strain energy density
	double l = m_E*m_v/((1+m_v)*(1-2*m_v));
	double m = 0.5*m_E/(1+m_v);
	
	double U = l*pow(e.tr(), 2)/2 + m*e.dotdot(e);
	
	return U;
}

//-----------------------------------------------------------------------------
//! calculate Cauchy stress
Matrix3ds FEElasticPlasticMaterial::Stress(FEMaterialPoint& pt)
{
	FEElasticPlasticMaterialPoint& et = *pt.ExtractData<FEElasticPlasticMaterialPoint>();
	
	// Get elastic strain
	Matrix3ds e = et.ElasticStrain();
	
	// Calculate elastic moduli
	double l = m_E*m_v/((1+m_v)*(1-2*m_v));
	double m = 0.5*m_E/(1+m_v);
	
	// Calculate trial stress
	Matrix3ds I(1.0);
	Matrix3ds s = I*(l*e.tr()) + e*(2*m);
	
	// Check for yielding
	double s_norm = sqrt(s.dotdot(s));
	double sy_scaled = m_sy; // In real implementation, this might be modified by hardening
	
	if (s_norm > sy_scaled) {
		// Plastic correction
		Matrix3ds s_dev = s - I*(s.tr()/3.0);
		double s_dev_norm = sqrt(s_dev.dotdot(s_dev));
		
		if (s_dev_norm > 0.0) {
			s = s_dev * (sy_scaled/s_dev_norm) + I*(s.tr()/3.0);
			// Mark as plastic
			et.m_bPlastic = true;
		}
	}
	
	// Store stress in material point
	et.m_s = s;
	
	return s;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness
tens4ds FEElasticPlasticMaterial::Tangent(FEMaterialPoint& pt)
{
	// For simplicity, we return the elastic tangent
	// A full implementation would modify this based on plastic state
	
	double l = m_E*m_v/((1+m_v)*(1-2*m_v));
	double m = 0.5*m_E/(1+m_v);
	
	tens4ds c;
	c.zero();
	Matrix3dd I(1.0);
	
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l_idx=0; l_idx<3; ++l_idx)
				{
					c(i,j,k,l_idx) = I(i,k)*I(j,l_idx)*l + (I(i,l_idx)*I(j,k) + I(i,k)*I(j,l_idx))*m;
				}
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate the 2nd P-K stress
Matrix3ds FEElasticPlasticMaterial::PK2Stress(FEMaterialPoint& pt, const Matrix3ds E)
{
	// Get the deformation gradient
	FEElasticPlasticMaterialPoint& et = *pt.ExtractData<FEElasticPlasticMaterialPoint>();
	Matrix3d F = et.m_F;
	double J = et.m_J;
	Matrix3d FiT = F.transinv();
	
	// Calculate Cauchy stress
	Matrix3ds s = Stress(pt);
	
	// Convert to 2nd P-K stress
	Matrix3d S = (FiT*(s*F.trans()))/J;
	Matrix3ds Spk2(S);
	
	return Spk2;
}

//-----------------------------------------------------------------------------
FEElasticPlasticStress::FEElasticPlasticStress() : FEDomainParameter("stress")
{

}

FEParamValue FEElasticPlasticStress::value(FEMaterialPoint& mp)
{
	FEElasticPlasticMaterialPoint& ep = *mp.ExtractData<FEElasticPlasticMaterialPoint>();
	return ep.m_s;
}