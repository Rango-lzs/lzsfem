#include "RgMaterialPointData.h"
#include "basicio/DumpStream.h"

namespace RgFem {
namespace Framework {

MaterialPointData::MaterialPointData(RgMaterialPointData* ppt)
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

void MaterialPointData::updateKinematicsFromF(const Matrix3d& F_new)
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

void MaterialPointData::pushForwardStress()
{
    // For finite strain, push forward the 2nd PK to Cauchy stress
    // sigma = 1/J * F * S * F^T
    sigma = (1.0 / J) * F * S * F.transpose();
}

void MaterialPointData::commit()
{
    ivar_committed = ivar_trial;
}

void MaterialPointData::revert()
{
    ivar_trial = ivar_committed;
}

void MaterialPointData::Init()
{
    RgMaterialPointData::Init();
    
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

void MaterialPointData::Update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::Update(timeInfo);
    
    // Update any time-dependent state variables here
    // For now, just copy current values to previous
    Fprev = F;
}

void MaterialPointData::Serialize(DumpStream& ar)
{
    RgMaterialPointData::Serialize(ar);
    
    // Serialize kinematic variables
    ar & F & Fprev & J & C & E & e;
    
    // Serialize stress measures
    ar & S & sigma;
    
    // Serialize internal variables
    ar & ivar_committed & ivar_trial;
}

} // namespace Framework

namespace LargeDef {

LargeDefMaterialPointData::LargeDefMaterialPointData(RgMaterialPointData* ppt)
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

void LargeDefMaterialPointData::updateKinematicsFromF(const Matrix3d& F_new)
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

void LargeDefMaterialPointData::pushForwardStress()
{
    // For finite strain, push forward the 2nd PK to Cauchy stress
    // sigma = 1/J * F * S * F^T
    sigma = (1.0 / J) * F * S * F.transpose();
}

void LargeDefMaterialPointData::commit()
{
    ivar_committed = ivar_trial;
}

void LargeDefMaterialPointData::revert()
{
    ivar_trial = ivar_committed;
}

void LargeDefMaterialPointData::Init()
{
    RgMaterialPointData::Init();
    
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

void LargeDefMaterialPointData::Update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::Update(timeInfo);
    
    // Update any time-dependent state variables here
    // For now, just copy current values to previous
    Fprev = F;
}

void LargeDefMaterialPointData::Serialize(DumpStream& ar)
{
    RgMaterialPointData::Serialize(ar);
    
    // Serialize kinematic variables
    ar & F & Fprev & J & C & E & e;
    
    // Serialize stress measures
    ar & S & sigma;
    
    // Serialize internal variables
    ar & ivar_committed & ivar_trial;
}

} // namespace LargeDef

namespace SmallDef {

SmallDefMaterialPointData::SmallDefMaterialPointData(RgMaterialPointData* ppt)
    : RgMaterialPointData(ppt)
{
    strain.zero();
    stress.zero();
    strain_prev.zero();
    
    elastic_strain.zero();
    plastic_strain.zero();
    eq_plastic_strain = 0.0;
    
    temperature = 0.0;
    damage = 0.0;
    
    state_vars.clear();
    state_vars_prev.clear();
}

void SmallDefMaterialPointData::Init()
{
    RgMaterialPointData::Init();
    
    // Initialize small deformation specific data
    strain.zero();          // Initialize strain to zero
    stress.zero();          // Initialize stress to zero
    strain_prev.zero();     // Initialize previous strain to zero
    
    elastic_strain.zero();      // Initialize elastic strain to zero
    plastic_strain.zero();      // Initialize plastic strain to zero
    eq_plastic_strain = 0.0;    // Initialize equivalent plastic strain to zero
    
    temperature = 0.0;      // Initialize temperature to reference value
    damage = 0.0;           // Initialize damage to zero (undamaged state)
    
    // Initialize state variables if needed
    state_vars.clear();
    state_vars_prev.clear();
}

void SmallDefMaterialPointData::Update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::Update(timeInfo);
    
    // Update any time-dependent state variables here
    // For small deformation, we might update strain_prev from strain
    // and update state variable history
    strain_prev = strain;
    state_vars_prev = state_vars;
}

void SmallDefMaterialPointData::Serialize(DumpStream& ar)
{
    RgMaterialPointData::Serialize(ar);
    
    // Serialize small deformation specific variables
    ar & strain & stress & strain_prev;
    ar & state_vars & state_vars_prev;
    ar & elastic_strain & plastic_strain & eq_plastic_strain;
    ar & temperature & damage;
}

} // namespace SmallDef

} // namespace RgFem