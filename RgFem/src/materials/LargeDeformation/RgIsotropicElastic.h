#pragma once
#include "femcore/RTTI/RTTIMacroDefine.h"
#include "materials/LargeDeformation/RgLargeDefMaterial.h"
#include "femcore/FEModelParam.h"

//此材料用于小应变，大转动场景
class FEIsotropicElastic : public LargeDef::RgLargeDefMaterial
{
    //DECLARE_META_CLASS(FEIsotropicElastic, LargeDef::RgLargeDefMaterial);

public:
    FEIsotropicElastic()
        : LargeDef::RgLargeDefMaterial()
    {
    }

    //! calculate stress at material point
    virtual Matrix3ds Stress(RgMaterialPoint& pt);

    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(RgMaterialPoint& pt);

    //! calculate strain energy density at material point
    virtual double StrainEnergyDensity(RgMaterialPoint& pt);

    //! calculate the 2nd Piola-Kirchhoff stress at material point
    Matrix3ds PK2Stress(RgMaterialPoint& pt, const Matrix3ds E);

    //! calculate material tangent stiffness at material point
    tens4dmm MaterialTangent(RgMaterialPoint& pt, const Matrix3ds E);

    // declare the parameter list
    //DECLARE_PARAM_LIST();

private:
    FEParamDouble m_E;  //!< Young's modulus
    FEParamDouble m_v;  //!< Poisson's ratio
};

