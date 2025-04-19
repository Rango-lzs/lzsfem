#pragma once
#include "materials/FEElasticMaterial.h"

//此材料用于小应变，大转动场景
class FEIsotropicElastic : public FEElasticMaterial
{
public:
    FEIsotropicElastic(FEModel* pfem)
        : FEElasticMaterial(pfem)
    {
    }

    //! calculate stress at material point
    virtual Matrix3ds Stress(FEMaterialPoint& pt) override;

    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(FEMaterialPoint& pt) override;

    //! calculate strain energy density at material point
    virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;

    //! calculate the 2nd Piola-Kirchhoff stress at material point
    Matrix3ds PK2Stress(FEMaterialPoint& pt, const Matrix3ds E) override;

    //! calculate material tangent stiffness at material point
    tens4dmm MaterialTangent(FEMaterialPoint& pt, const Matrix3ds E) override;

    // declare the parameter list
    DECLARE_PARAM_LIST();

private:
    FEParamDouble m_E;  //!< Young's modulus
    FEParamDouble m_v;  //!< Poisson's ratio
};

