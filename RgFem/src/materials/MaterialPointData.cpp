#include "framework/MaterialPointData.h"
#include "datastructure/mathalg.h"

namespace RgFem {
namespace Framework {

MaterialPointData::MaterialPointData(FEMaterialPointData* ppt)
    : FEMaterialPointData(ppt)
    , J(1.0)
{
    // Initialize identity tensors
    F.identity();
    Fprev.identity();
    C.identity();
    E.zero();
    e.zero();
    S.zero();
    sigma.zero();
}

void MaterialPointData::updateKinematicsFromF(const Matrix3d& F_new)
{
    F = F_new;
    J = F.det();
    
    // Compute right Cauchy-Green tensor: C = F^T * F
    C = F.transpose() * F;
    
    // Compute Green-Lagrange strain: E = 0.5 * (C - I)
    E = C;
    E[0][0] -= 1.0;
    E[1][1] -= 1.0;
    E[2][2] -= 1.0;
    E *= 0.5;
    
    // For small strain, we might also store the infinitesimal strain
    // This would be computed differently based on the analysis type
}

void MaterialPointData::pushForwardStress()
{
    // Push forward PK2 stress to Cauchy stress:
    // sigma = (1/J) * F * S * F^T
    if (J > 0.0) {
        Matrix3d FT = F.transpose();
        sigma = (1.0 / J) * (F * S * FT);
    } else {
        sigma.zero();
    }
}

void MaterialPointData::commit()
{
    ivar_committed = ivar_trial;
    Fprev = F;
}

void MaterialPointData::revert()
{
    ivar_trial = ivar_committed;
    F = Fprev;
    updateKinematicsFromF(F);
}

void MaterialPointData::init()
{
    FEMaterialPointData::init();
    // Initialize to identity
    F.identity();
    Fprev.identity();
    J = 1.0;
    C.identity();
    E.zero();
    e.zero();
    S.zero();
    sigma.zero();
}

void MaterialPointData::update(const FETimeInfo& timeInfo)
{
    FEMaterialPointData::update(timeInfo);
    // Update any time-dependent variables if needed
}

void MaterialPointData::serialize(DumpStream& ar)
{
    FEMaterialPointData::serialize(ar);
    
    if (ar.IsSaving()) {
        ar << J;
        ar << F << Fprev << C << E << e;
        ar << S << sigma;
        ar << ivar_committed << ivar_trial;
    } else {
        ar >> J;
        ar >> F >> Fprev >> C >> E >> e;
        ar >> S >> sigma;
        ar >> ivar_committed >> ivar_trial;
    }
}

} // namespace Framework
} // namespace RgFem