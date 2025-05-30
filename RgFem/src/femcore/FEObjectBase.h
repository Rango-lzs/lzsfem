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
#include "femcore/FEProperty.h"

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
    //! Add a property
    //! Call this in the constructor of derived classes to
    //! build the property list
    
    void AddProperty( FEProperty* pp, const char* sz, unsigned int flags);

    //-----------------------------------------------------------------------------
    void RemoveProperty(int i);

    //-----------------------------------------------------------------------------
    void ClearProperties();

    //-----------------------------------------------------------------------------
    int Properties();

    //-----------------------------------------------------------------------------
    int FindPropertyIndex(const char* sz);

    //-----------------------------------------------------------------------------
    FEProperty* FindProperty(const char* sz, bool searchChildren = false);
    
    //-----------------------------------------------------------------------------
    FEProperty* GetProperty(int n);


    //-----------------------------------------------------------------------------
    bool SetProperty(int i, FEObjectBase* pb);

    //-----------------------------------------------------------------------------
    //! Set a property via name
    bool SetProperty(const char* sz, FEObjectBase* pb);

    template <class T>
    T* ExtractProperty(bool extractSelf = true);

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
    std::vector<FEProperty*>	m_Prop;		//!< list of properties
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


// include template property definitions
#include "FEPropertyT.h"

template <class T>
FEProperty* AddClassProperty(FEObjectBase* pc, T* pp, const char* sz)
{
    FEFixedPropertyT<T>* prop = new FEFixedPropertyT<T>(pp);
    prop->SetDefaultType(sz);
    pc->AddProperty(prop, sz, FEProperty::Fixed);
    return prop;
}

template <class T>
FEProperty* AddClassProperty(FEObjectBase* pc, T** pp, const char* sz, unsigned int flags = FEProperty::Required)
{
    FEPropertyT<T>* prop = new FEPropertyT<T>(pp);
    if (prop->GetSuperClassID() == FECLASS_ID)
        prop->SetDefaultType(sz);
    pc->AddProperty(prop, sz, flags);
    return prop;
}

template <class T>
FEProperty* AddClassProperty(FEObjectBase* pc, std::vector<T*>* pp, const char* sz,
                             unsigned int flags = FEProperty::Required)
{
    FEVecPropertyT<T>* prop = new FEVecPropertyT<T>(pp);
    if (prop->GetSuperClassID() == FECLASS_ID)
        prop->SetDefaultType(sz);
    pc->AddProperty(prop, sz, flags);
    return prop;
}

#define ADD_PROPERTY(theProp, ...) AddClassProperty(this, &theProp, __VA_ARGS__)

template <class T>
T* FEObjectBase::ExtractProperty(bool extractSelf)
{
    if (extractSelf)
    {
        if (dynamic_cast<T*>(this))
            return dynamic_cast<T*>(this);
    }

    int NC = Properties();
    for (int i = 0; i < NC; i++)
    {
        FEObjectBase* pci = GetProperty(i);
        if (pci)
        {
            T* pc = pci->ExtractProperty<T>();
            if (pc)
                return pc;
        }
    }
    return nullptr;
}
