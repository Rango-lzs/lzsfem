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
#include <functional>
#include <string>

// �Ƿ���Ҫһ��RxStore���������е�RxClass��Ϣ��

class MetaObject;
using ObjectConstructor = std::function<MetaObject*()>;

class MetaClass
{
public:
    // ��Ҫͨ�����캯����������MetaClass���ܹ�ʵ����MetaObject�Ĺ��캯��
    MetaClass(std::string name, const MetaClass* pParent, ObjectConstructor cons);
    ~MetaClass() = default;

    const std::string& name() const;
    const MetaClass* parent() const;
    bool isKindOf(const MetaClass* pParent) const;

    MetaObject* create();

private:
    std::string m_name;
    const MetaClass* mp_parent;
    ObjectConstructor m_constructor;
};

template <class T>
class ConcretMeta : public MetaClass
{
public:
    static const MetaClass* instance()
    {
        static ConcretMeta s_meta;
        MetaClassStore::instance()->insert(&s_meta);
        return &s_meta;
    }

private:
    ConcretMeta()
        : MetaClass("name", ConcretMeta<T::BaseClass>::instance(), []() ->  MetaObject* { return new T(); })
    {
    }
};

#endif