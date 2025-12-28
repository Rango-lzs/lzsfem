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
    : RgMaterialPointData(ppt), gradU(0), strain(0), strain_prev(0), 
      stress(0), stress_prev(0), C(0) 
{
}

void SmallDefMaterialPointData::updateKinematicsFromGradU(const Matrix3d& gradU_new)
{
    gradU = gradU_new;
    
    // Calculate infinitesimal strain: ε = 1/2(∇u + ∇u^T)
    Matrix3d temp_strain = 0.5 * (gradU + gradU.transpose());
    // Convert to symmetric tensor
    strain = temp_strain.sym();
}

void SmallDefMaterialPointData::commit()
{
    strain_prev = strain;
    stress_prev = stress;
    ivar_committed = ivar_trial;
}

void SmallDefMaterialPointData::revert()
{
    strain = strain_prev;
    stress = stress_prev;
    ivar_trial = ivar_committed;
}

void SmallDefMaterialPointData::init()
{
    // Initialize to zero or identity as appropriate
    Matrix3d zero_matrix(0, 0, 0, 0, 0, 0, 0, 0, 0);
    Matrix3ds zero_sym_matrix(0, 0, 0, 0, 0, 0);  // xx, yy, zz, xy, yz, xz
    
    gradU = zero_matrix;
    strain = zero_sym_matrix;
    strain_prev = zero_sym_matrix;
    stress = zero_sym_matrix;
    stress_prev = zero_sym_matrix;
    C = tens4d(0);  // Initialize fourth-order tensor to zero
    
    ivar_committed.clear();
    ivar_trial.clear();
    
    RgMaterialPointData::init();
}

void SmallDefMaterialPointData::update(const FETimeInfo& timeInfo)
{
    RgMaterialPointData::update(timeInfo);
}

void SmallDefMaterialPointData::serialize(DumpStream& ar)
{
    ar & gradU & strain & strain_prev & stress & stress_prev & ivar_committed & ivar_trial & C;
    RgMaterialPointData::serialize(ar);
}

} // namespace SmallDef

} // namespace RgFem