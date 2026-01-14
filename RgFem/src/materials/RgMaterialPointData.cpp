#include "RgMaterialPointData.h"
#include "basicio/DumpStream.h"

namespace LargeDef {

LargeDefRgMaterialPointData::LargeDefRgMaterialPointData(RgMaterialPointData* ppt)
    : RgMaterialPointData(ppt)
{
    F.zero();
    Fprev.unit();
    J = 1.0;
    C.unit();
    E.zero();
    e.zero();
    
    S.zero();
    sigma.zero();
    
    ivar_committed.clear();
    ivar_trial.clear();
}

void LargeDefRgMaterialPointData::updateKinematicsFromF(const Matrix3d& F_new)
{
    F = F_new;
    J = F.det();
    C = F.transpose() * F;  // Right Cauchy-Green: C = F^T * F
    E = 0.5 * (C - Matrix3d::identity());  // Green-Lagrange strain: E = 0.5(C - I)
    
    // Small strain (for small deformation assumption)
    // This is typically the symmetric part of the displacement gradient
    Matrix3d F_sym = 0.5 * (F + F.transpose());
    e = 0.5 * (Matrix3d::identity() - F_sym.inverse().transpose());  // Almansi strain
}

void LargeDefRgMaterialPointData::pushForwardStress()
{
    // For finite strain, push forward the 2nd PK to Cauchy stress
    // sigma = 1/J * F * S * F^T
    sigma = (1.0 / J) * F * S * F.transpose();
}

void LargeDefRgMaterialPointData::commit()
{
    ivar_committed = ivar_trial;
}

void LargeDefRgMaterialPointData::revert()
{
    ivar_trial = ivar_committed;
}

void LargeDefRgMaterialPointData::init()
{
    RgMaterialPointData::init();
    
    F.unit();           // Initialize to identity
    Fprev.unit();
    J = 1.0;           // Initialize to incompressible (det F = 1)
    C.unit();
    E.zero();
    e.zero();
    
    S.zero();
    sigma.zero();
    
    // Initialize internal variables if needed
    ivar_committed.clear();
    ivar_trial.clear();
}

void LargeDefRgMaterialPointData::update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::update(timeInfo);
    
    // Update any time-dependent state variables here
    // For now, just copy current values to previous
    Fprev = F;
}

void LargeDefRgMaterialPointData::serialize(DumpStream& ar)
{
    RgMaterialPointData::serialize(ar);
    
    // Serialize kinematic variables
    ar & F & Fprev & J & C & E & e;
    
    // Serialize stress measures
    ar & S & sigma;
    
    // Serialize internal variables
    ar & ivar_committed & ivar_trial;
}

} // namespace LargeDef

namespace SmallDef {

SmallDefRgMaterialPointData::SmallDefRgMaterialPointData(RgMaterialPointData* ppt)
    : RgMaterialPointData(ppt)
{
    gradU.zero();
    strain.zero();
    strain_prev.zero();
    stress.zero();
    stress_prev.zero();
    C.zero();   
}

void SmallDefRgMaterialPointData::updateKinematicsFromGradU(const Matrix3d& gradU_new)
{
    gradU = gradU_new;
    
    // Calculate infinitesimal strain: ε = 1/2(∇u + ∇u^T)
    strain = 0.5 * (gradU + gradU.transpose());
}

void SmallDefRgMaterialPointData::commit()
{
    strain_prev = strain;
    stress_prev = stress;
}

void SmallDefRgMaterialPointData::revert()
{
    strain = strain_prev;
    stress = stress_prev;
}

void SmallDefRgMaterialPointData::init()
{
    RgMaterialPointData::init();
    
    // Initialize small deformation specific data
    gradU.zero();           // Initialize displacement gradient to zero
    strain.zero();          // Initialize strain to zero
    strain_prev.zero();     // Initialize previous strain to zero
    stress.zero();          // Initialize stress to zero
    stress_prev.zero();     // Initialize previous stress to zero
    C.zero();               // Initialize elasticity tensor to zero 
}

void SmallDefRgMaterialPointData::update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::update(timeInfo);
    
    // Update any time-dependent state variables here
    // For small deformation, we update the previous values during commit
    // This is handled in commit() method
}

void SmallDefRgMaterialPointData::serialize(DumpStream& ar)
{
    RgMaterialPointData::serialize(ar);
    
    // Serialize small deformation specific variables
    ar & gradU & strain & strain_prev & stress & stress_prev;
}

} // namespace SmallDef

