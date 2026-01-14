#include "materials/LargeDeformation/RgIsotropicElastic.h"
#include "femcore/units.h"
#include "femcore/FEParamValidator.h"

DEFINE_META_CLASS(FEIsotropicElastic, FEElasticMaterial, "iso-elastic");

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAM_DEFINE(FEIsotropicElastic, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v")->setLongName("Poisson's ratio");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
Matrix3ds FEIsotropicElastic::Stress(RgMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	Matrix3d &F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = Ji*(v*E/((1+v)*(1-2*v)));
	double mu  = Ji*(0.5*E/(1+v));

	// calculate left Cauchy-Green tensor (ie. b-Matrix)
	Matrix3ds b = pt.LeftCauchyGreen();

	// calculate trace of Green-Lagrance strain tensor
	double trE = 0.5*(b.tr()-3);

	// calculate square of b-Matrix
	// (we commented out the Matrix components we do not need)
	Matrix3ds b2 = b.sqr();

	// calculate stress
	Matrix3ds s = b*(lam*trE - mu) + b2*mu;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEIsotropicElastic::Tangent(RgMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double E =  m_E(mp);
    double v = 0.3; //m_v(mp);

	// deformation gradient
	Matrix3d& F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	// lame parameters
	double lam = Ji*(v*E/((1+v)*(1-2*v)));
	double mu  = Ji*(0.5*E/(1+v));

	// left cauchy-green Matrix (i.e. the 'b' Matrix)
	Matrix3ds b = pt.LeftCauchyGreen();

	return dyad1s(b)*lam + dyad4s(b)*(2.0*mu);
}

//-----------------------------------------------------------------------------
double FEIsotropicElastic::StrainEnergyDensity(RgMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double mE = m_E(mp);
	double mv = m_v(mp);

    Matrix3ds E = (pt.RightCauchyGreen() - Matrix3dd(1))/2;
    
    double lam = mv*mE/((1+mv)*(1-2*mv));
	double mu  = 0.5*mE/(1+mv);

    double trE = E.tr();
    double Enorm = E.norm();
    
    double sed = lam*trE*trE/2 + mu*Enorm*Enorm;
    
	return sed;
}

//-----------------------------------------------------------------------------
Matrix3ds FEIsotropicElastic::PK2Stress(RgMaterialPoint& pt, const Matrix3ds E)
{
	double mE = m_E(pt);
	double mv = m_v(pt);

    // lame parameters
    double lam = (mv*mE/((1+mv)*(1-2*mv)));
    double mu  = (0.5*mE/(1+mv));
    
    // Identity
    Matrix3dd I(1);
    
    // calculate stress
    Matrix3ds S = I*(E.tr()*lam) + E*(2*mu);
    
    return S;
}

//-----------------------------------------------------------------------------
tens4dmm FEIsotropicElastic::MaterialTangent(RgMaterialPoint& pt, const Matrix3ds E)
{
	double mE = m_E(pt);
	double mv = m_v(pt);

    // lame parameters
    double lam = (mv*mE/((1+mv)*(1-2*mv)));
    double mu  = (0.5*mE/(1+mv));
    
    // Identity
    Matrix3dd I(1);
    
    tens4dmm c = dyad1s(I)*lam + dyad4s(I)*(2*mu);
    
    return c;
}
