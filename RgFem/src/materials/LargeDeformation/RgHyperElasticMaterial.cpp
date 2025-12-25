#include "FEElasticMaterial.h"
#include "femcore/FEModel.h"

DEFINE_META_CLASS(FEElasticMaterial, FESolidMaterial, "");

BEGIN_PARAM_DEFINE(FEElasticMaterial, FESolidMaterial)
END_PARAM_DEFINE();

FEElasticMaterial::FEElasticMaterial() : FESolidMaterial()
{ 
	m_density = 1;
	//AddDomainParameter(new FEElasticStress());
}

//-----------------------------------------------------------------------------
FEElasticMaterial::~FEElasticMaterial()
{ 
	
}

//-----------------------------------------------------------------------------
RgMaterialPointData* FEElasticMaterial::CreateMaterialPointData()
{ 
	return new FEElasticMaterialPoint;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
Matrix3ds FEElasticMaterial::SecantStress(FEMaterialPoint& mp, bool PK2)
{
	// extract the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	Matrix3d F = pt.m_F;
	double J = pt.m_J;
	Matrix3ds E = pt.Strain();
	Matrix3dd I(1);
	Matrix3d FiT = F.transinv();

	// calculate the 2nd P-K stress at the current deformation gradient
	double W = StrainEnergyDensity(mp);

	// create deformation gradient increment
	double eps = 1e-9;
	Vector3d e[3];
	e[0] = Vector3d(1, 0, 0); e[1] = Vector3d(0, 1, 0); e[2] = Vector3d(0, 0, 1);
	Matrix3ds S(0.0);
	for (int k = 0; k < 3; ++k) {
		for (int l = k; l < 3; ++l) {
			// evaluate incremental stress
			Matrix3d dF = FiT * ((e[k] & e[l]))*(eps*0.5);
			Matrix3d F1 = F + dF;
			pt.m_F = F1;
			pt.m_J = pt.m_F.det();

			double dW = StrainEnergyDensity(mp) - W;

			// evaluate the secant modulus
			S(k, l) = 2.0 * dW / eps;
		}
	}

    // restore values
    pt.m_F = F;
    pt.m_J = J;
    
    if (PK2) return S;
    else {
        // push from material to spatial frame
        Matrix3ds s = pt.push_forward(S);
        
        // return secant stress
        return s;
    }
}

//-----------------------------------------------------------------------------
//! return the strain energy density
double FEElasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt) { return 0; }

//-----------------------------------------------------------------------------
FEElasticStress::FEElasticStress() : FEDomainParameter("stress")
{

}

FEParamValue FEElasticStress::value(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	return ep.m_s;
}