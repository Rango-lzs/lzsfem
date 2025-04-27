/*****************************************************************
 * \file   meta_object.h
 * \brief
 *
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once
#ifndef META_OBJECT_H
#define META_OBJECT_H

#include "femcore/fem_export.h"

#include <string>

class MetaClass;
class FEModel;

class FEM_EXPORT MetaObject
{
public:
    using BaseClass = void;  // 这个需要单独定义，否则其子类MetaClass构造函数找不到BaseClass
    virtual ~MetaObject() = 0;

    virtual const MetaClass* meta() const = 0;  // called by instance
    static const MetaClass* staic_meta();       // called by class
    static MetaObject* meta_cast(MetaObject* pOther);
    static std::string class_name();
    bool isKindOf(const MetaClass* pMeta) const;

public:
    MetaObject() = default;
};

// 遍历
template <class T>
T* RANGO_NEW(const std::string& aliasName = "")
{
    const MetaClass* pTarget = nullptr;
    T* pNewObj = nullptr;
    if (aliasName.empty() || aliasName == T::alias_name())
    {
        pTarget = T::static_meta();
    }
    else
    {
        for (auto pMeta : pTarget->m_childs)
        {
            if (pMeta->m_alias_name == aliasName)
            {
                pTarget = pMeta;
            }
        }
    }
    pNewObj = static_cast<T*>(pMeta->create());
    return pNewObj;
}

template <class T>
T* RANGO_NEW(FEModel* pModel, const std::string& aliasName)
{
    const MetaClass* pTarget = nullptr;
    T* pNewObj = nullptr;
    if (aliasName.empty() || aliasName == T::alias_name())
    {
        pTarget = T::static_meta();
    }
    else
    {
        for (auto pMeta : pTarget->m_childs)
        {
            if (pMeta->m_alias_name == aliasName)
            {
                pTarget = pMeta;
            }
        }
    }
    pNewObj = static_cast<T*>(pMeta->create());
    pNewObj->SetModel(pModel);
    return pNewObj;
}

#endif