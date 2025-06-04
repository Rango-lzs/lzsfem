/*********************************************************************
 * \file   FESolidMaterial.cpp
 * \brief  
 * 
 * \author Leizs
 * \date   April 2025
 *********************************************************************/

#include "materials/FESolidMaterial.h"
#include "materials/FEElasticMaterial.h"
#include "femcore/units.h"

// Material parameters for FEElasticMaterial
BEGIN_PARAM_DEFINE(FESolidMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density")->setUnits(UNIT_DENSITY)->MakeTopLevel(true);
END_PARAM_DEFINE();

FESolidMaterial::FESolidMaterial() : FEMaterial()
{
    m_density = 1.0;
}

//! set the material density
void FESolidMaterial::SetDensity(double d)
{ 
	m_density = d;
}

//! evaluate density
double FESolidMaterial::Density(FEMaterialPoint& pt)
{
	return m_density(pt);
}

tens4dmm FESolidMaterial::SolidTangent(FEMaterialPoint& mp)
{
	return (UseSecantTangent() ? SecantTangent(mp) : Tangent(mp));
}

//-----------------------------------------------------------------------------
Matrix3ds FESolidMaterial::SecantStress(FEMaterialPoint& pt, bool PK2)
{
    assert(false);
    return Matrix3ds(0.0);
}

//-----------------------------------------------------------------------------
//! calculate the 2nd Piola-Kirchhoff stress at material point, using prescribed Lagrange strain
//! needed for EAS analyses where the compatible strain (calculated from displacements) is enhanced
Matrix3ds FESolidMaterial::PK2Stress(FEMaterialPoint& mp, const Matrix3ds E)
{
    // Evaluate right Cauchy-Green tensor from E
    Matrix3ds C = Matrix3dd(1) + E*2;
    
    // Evaluate right stretch tensor U from C
    Vector3d v[3];
    double lam[3];
    C.eigen2(lam, v);
    lam[0] = sqrt(lam[0]); lam[1] = sqrt(lam[1]); lam[2] = sqrt(lam[2]);
    Matrix3ds U = dyad(v[0])*lam[0] + dyad(v[1])*lam[1] + dyad(v[2])*lam[2];
    double J = lam[0]*lam[1]*lam[2];
    
    // temporarily replace F in material point with U
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    Matrix3d Fsafe = pt.m_F;
    double Jsafe = pt.m_J;
    pt.m_F = U;
    pt.m_J = J;
    
    // Evaluate Cauchy stress
    Matrix3ds s = Stress(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert Cauchy stress to 2nd P-K stress
    Matrix3ds Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    Matrix3ds S = (Ui*s*Ui).sym()*J;
    
    return S;
}

//-----------------------------------------------------------------------------
//! calculate material tangent stiffness at material point, using prescribed Lagrange strain
//! needed for EAS analyses where the compatible strain (calculated from displacements) is enhanced
tens4dmm FESolidMaterial::MaterialTangent(FEMaterialPoint& mp, const Matrix3ds E)
{
    // Evaluate right Cauchy-Green tensor from E  E= 1/2(C-I)
    Matrix3ds C = Matrix3dd(1) + E*2;
    
    // Evaluate right stretch tensor U from C
    Vector3d v[3];
    double lam[3];
    C.eigen2(lam, v);
    lam[0] = sqrt(lam[0]); lam[1] = sqrt(lam[1]); lam[2] = sqrt(lam[2]);
    Matrix3d U = dyad(v[0])*lam[0] + dyad(v[1])*lam[1] + dyad(v[2])*lam[2];
    double J = lam[0]*lam[1]*lam[2];
    
    // temporarily replace F in material point with U
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    Matrix3d Fsafe = pt.m_F;
    double Jsafe = pt.m_J;
    pt.m_F = U;
    pt.m_J = J;
    
    // Evaluate Cauchy stress
    tens4dmm c = SolidTangent(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert spatial tangent to material tangent
    Matrix3d Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    tens4dmm Cm = c.pp(Ui)*J;
    
    return Cm;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
tens4dmm FESolidMaterial::SecantTangent(FEMaterialPoint& mp, bool mat)
{
    // extract the deformation gradient
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    Matrix3d F = pt.m_F;
    double J = pt.m_J;
    Matrix3ds E = pt.Strain();
    Matrix3dd I(1);

    // calculate the 2nd P-K stress at the current deformation gradient
    Matrix3ds S = PK2Stress(mp,E);
    
    // create deformation gradient increment
    double eps = 1e-9;
    Vector3d e[3];
    e[0] = Vector3d(1,0,0); e[1] = Vector3d(0,1,0); e[2] = Vector3d(0,0,1);
    tens4dmm C;
    for (int k=0; k<3; ++k) {
        for (int l=k; l<3; ++l) {
            // evaluate incremental stress
            Matrix3ds dE = dyads(e[k], e[l])*(eps/2);
            Matrix3ds dS = (PK2Stress(mp,E+dE) - S)/eps;
            
            // evaluate the secant modulus
            C(0,0,k,l) = C(0,0,l,k) = dS.xx();
            C(1,1,k,l) = C(1,1,l,k) = dS.yy();
            C(2,2,k,l) = C(2,2,l,k) = dS.zz();
            C(0,1,k,l) = C(1,0,k,l) = C(0,1,l,k) = C(1,0,l,k) = dS.xy();
            C(1,2,k,l) = C(2,1,k,l) = C(1,2,l,k) = C(2,1,l,k) = dS.yz();
            C(2,0,k,l) = C(0,2,k,l) = C(2,0,l,k) = C(0,2,l,k) = dS.xz();
        }
    }
    
    if (mat) return C;
    else {
        
        // push from material to spatial frame
        tens4dmm c = C.pp(F)/J;
        
        // return secant tangent
        return c;
    }
}
