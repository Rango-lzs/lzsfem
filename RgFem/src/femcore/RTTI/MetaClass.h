/*****************************************************************
 * \file   MetaClass.h
 * \brief  
 * 
 * \author 11914
 * \date   March 2025
 *********************************************************************/
#pragma once

#ifndef META_CLASS_H
#define META_CLASS_H

#include "femcore/RTTI/MetaClassStore.h"
#include "femcore/MetaObject.h"
#include <functional>
#include <string>
#include <list>
#include <type_traits>

// 是否需要一个RxStore来管理所有的RxClass信息？

class MetaObject;
using ObjectConstructor = std::function<MetaObject*()>;

class MetaClass
{
public:
    // 需要通过构造函数传入基类的MetaClass和能够实例化MetaObject的构造函数
    MetaClass(const std::string& className, const std::string& aliasName, MetaClass* pParent, ObjectConstructor cons);
    ~MetaClass() = default;

    const std::string& name() const;
    const MetaClass* parent() const;
    bool isKindOf(const MetaClass* pParent) const;
    MetaObject* create() const;
    
private:
    std::string m_class_name;
    std::string m_alias_name;
    const MetaClass* mp_parent;
    std::list<MetaClass*> m_childs;
    ObjectConstructor m_constructor;
};

template<typename T, typename Y = void>
class ConcreteMeta;

template<typename T>
class ConcreteMeta<typename T, typename std::enable_if_t<!std::is_abstract_v<T>, void> > : public MetaClass
{
public:
    static MetaClass* instance()
    {
        static ConcreteMeta s_meta;
        MetaClassStore::instance()->insert(&s_meta);
        return &s_meta;
    }

private:
    ConcreteMeta()
        : MetaClass(T::class_name(),T::alias_name(), ConcreteMeta<T::BaseClass>::instance(), []() -> MetaObject* { return new T(); })
    {
    }
};

template<typename T>
class ConcreteMeta<typename T, typename std::enable_if_t<std::is_abstract<T>::value, void>> : public MetaClass
{
public:
    static MetaClass* instance()
    {
        static ConcreteMeta s_meta;
        MetaClassStore::instance()->insert(&s_meta);
        return &s_meta;
    }

private:
    ConcreteMeta()
        : MetaClass(T::class_name(), T::alias_name(), ConcreteMeta<T::BaseClass>::instance(),
                    []() -> MetaObject* { return nullptr; })
    {
    }
};

template <>
class ConcreteMeta<MetaObject> : public MetaClass
{
public:
    static MetaClass* instance()
    {
        static ConcreteMeta s_meta;
        MetaClassStore::instance()->insert(&s_meta);
        return &s_meta;
    }

private:
    ConcreteMeta()
        : MetaClass(MetaObject::class_name(),"meta", nullptr, []() -> MetaObject* { return nullptr; })
    {
    }
};

#endif