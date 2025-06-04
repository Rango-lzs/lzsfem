#pragma once
#include "femcore/Domain/FEDomainList.h"
#include "femcore/Domain/FEDomainParameter.h"
#include "femcore/FEModelComponent.h"
#include "femcore/FEModelParam.h"
#include "datastructure/Matrix3d.h"

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEDomain;
class DumpStream;
class FEMaterialPoint;
class FEMaterialPointData;
class FEMat3dValuator;

//-----------------------------------------------------------------------------
class FEM_EXPORT FEMaterialBase : public FEModelComponent
{
    DECLARE_META_CLASS(FEMaterialBase, FEModelComponent);
public:
    FEMaterialBase();

    //! returns a pointer to a new material point object
    virtual FEMaterialPointData* CreateMaterialPointData() = 0;

    //! Update specialized material points at each iteration
    virtual void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) = 0;

    // evaluate local coordinate system at material point
    virtual Matrix3d GetLocalCS(const FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types
//! From this class all other material classes are derived.

class FEM_EXPORT FEMaterial : public FEMaterialBase
{
    DECLARE_META_CLASS(FEMaterial, FEMaterialBase);

public:

    FEMaterial();
    virtual ~FEMaterial();

    //! performs initialization
    bool Init() override;

    // evaluate local coordinate system at material point
    Matrix3d GetLocalCS(const FEMaterialPoint& mp) override;

    // set the (local) material axis valuator
    void SetMaterialAxis(FEMat3dValuator* val);

protected:
    FEMat3dValuator* m_Q;  //!< local material coordinate system

public:
    //! Assign a domain to this material
    void AddDomain(FEDomain* dom);

    //! get the domaint list
    FEDomainList& GetDomainList()
    {
        return m_domList;
    }

private:
    FEDomainList m_domList;                   //!< list of domains that use this material
    DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
// Material properties are classes that can only be defined as properties of other materials
class FEM_EXPORT FEMaterialProperty : public FEMaterialBase
{
public:
    FEMaterialProperty();

    // evaluate local coordinate system at material point
    Matrix3d GetLocalCS(const FEMaterialPoint& mp) override;

    DECLARE_PARAM_LIST();
};
