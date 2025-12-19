#include "materials/LargeDeformation/RgNeoHookean.h"
#include "datastructure/mathalg.h"
#include "materials/MaterialPointData.h"

namespace RgFem
{
    namespace LargeDef
    {

        NeoHookean::NeoHookean(double mu, double lambda)
            : m_mu(mu)
            , m_lambda(lambda)
        {
        }

        FEMaterialPointData* NeoHookean::createMaterialPointData() const
        {
            return new MaterialPointData();
        }

        void NeoHookean::computeConstitutive(MaterialPointData* mp, Matrix& D)
        {
            // For Neo-Hookean material, the constitutive relationship is nonlinear
            // So we would typically return the tangent matrix at the current state
            // This is a simplified implementation - a full implementation would
            // involve more detailed computation based on the current deformation state

            // For demonstration purposes, we're returning a zero matrix
            // A real implementation would compute the spatial tangent stiffness
            D.resize(6, 6);
            D.zero();

            // Store computed PK2 stress in material point data
            // mp->S = computePK2FromF(mp->F, mp->J);
        }

        void NeoHookean::commitState(MaterialPointData* mp)
        {
            mp->commit();
        }

        void NeoHookean::revertState(MaterialPointData* mp)
        {
            mp->revert();
        }

        std::string NeoHookean::getName() const
        {
            return "NeoHookean";
        }

        Matrix3d NeoHookean::computePK2FromF(const Matrix3d& F, double J) const
        {
            // Compute right Cauchy-Green tensor
            Matrix3d C = F.transpose() * F;

            // Compute inverse of C
            Matrix3d C_inv = C.inverse();

            // Compute PK2 stress for Neo-Hookean material:
            // S = mu * (I - C^-1) + lambda * (J - 1) * J * C^-1
            Matrix3d I;
            I.identity();

            Matrix3d S = m_mu * (I - C_inv);
            S += m_lambda * (J - 1.0) * J * C_inv;

            return S;
        }

    }  // namespace LargeDef
}  // namespace RgFem