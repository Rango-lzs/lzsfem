/*********************************************************************
 * \file   FEObjectBase.h
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/
#pragma once

#include "femcore/FEParamObject.h"
#include "femcore/RTTI/MetaClass.h"
#include <string>

//-----------------------------------------------------------------------------
class FECoreFactory;
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all
//! classes that can be registered with the framework.
class FEM_EXPORT FEObjectBase : public FEParamObject
{
    DECLARE_META_CLASS(FEObjectBase, FEParamObject);
public:
    //! constructor
    FEObjectBase();
    explicit FEObjectBase(FEModel* pModel);

    //! destructor
    virtual ~FEObjectBase();

    //! Initialization
    virtual bool Init();

    void SetName(const std::string& name);
    const std::string& GetName() const;

    virtual std::string GetTypeStr();

    void SetID(int id);
    int GetID() const;

public:
    //! number of parameters
    int Parameters() const;
    //! return a parameter
    virtual FEParam* FindParameter(const std::string& s) override;
    //! return the property (or this) that owns a parameter
    FEObjectBase* FindParameterOwner(void* pd);
    //! validates all properties and parameters
    bool Validate() override;
    //! call this after the parameters are changed
    virtual bool UpdateParams();

public:
    //! Get the FE model
    FEModel* GetFEModel() const;
    //! set the FEModel of this class (use with caution!)
    void SetFEModel(FEModel* fem);

    static void SaveClass(DumpStream& ar, FEObjectBase* p);
    static FEObjectBase* LoadClass(DumpStream& ar, FEObjectBase* p);
    //! data serialization
    void Serialize(DumpStream& ar) override;

private:
    std::string m_name;  //!< user defined name of component
    FEModel* m_fem;      //!< the model this class belongs to
    int m_nID;           //!< component ID
};


template <class T>
T* RANGO_NEW(const std::string& aliasName = "")
{
    const MetaClass* pThis = T::static_meta();
    const MetaClass* pTarget = nullptr;
    T* pNewObj = nullptr;
    if (aliasName.empty() || aliasName == T::alias_name())
    {
        pTarget = pThis;
    }
    else
    {
        for (auto pMeta : pThis->childs())
        {
            if (pMeta->alias_name() == aliasName)
            {
                pTarget = pMeta;
            }
        }
    }
    pNewObj = static_cast<T*>(pTarget->create());
    return pNewObj;
}

template <class T>
T* RANGO_NEW(FEModel* pModel, const std::string& aliasName)
{
    const MetaClass* pThis = T::static_meta();
    const MetaClass* pTarget = nullptr;
    T* pNewObj = nullptr;
    if (aliasName.empty() || aliasName == T::alias_name())
    {
        pTarget = pThis;
    }
    else
    {
        for (auto pMeta : pThis->childs())
        {
            if (pMeta->alias_name() == aliasName)
            {
                pTarget = pMeta;
            }
        }
    }
    pNewObj = static_cast<T*>(pTarget->create());
    pNewObj->SetFEModel(pModel);
    return pNewObj;
}