#pragma once

#ifndef META_CLASS_H
#define META_CLASS_H

#include <functional>
#include <string>
#include <type_traits>

class MetaObject;

//---------------------- ���ͱ������� ----------------------
using ObjectConstructor = std::function<MetaObject*()>;

//---------------------- �궨�壨���û�������----------------------
#define DECLARE_META_CLASS(ClassName, ParentClass)                                                                     \
public:                                                                                                                \
    using MetaBase = ParentClass;                                                                                      \
    static constexpr const char* StaticClassName()                                                                     \
    {                                                                                                                  \
        return #ClassName;                                                                                             \
    }                                                                                                                  \
    virtual const MetaClass* getMetaClass() const override                                                             \
    {                                                                                                                  \
        return ConcretMeta<ClassName>::instance();                                                                     \
    }                                                                                                                  \
                                                                                                                       \
private:

//---------------------- MetaClass ���Ķ��� ----------------------
class MetaClass
{
public:
    MetaClass(const std::string& name, const MetaClass* pParent, ObjectConstructor cons)
        : m_name(name)
        , mp_parent(pParent)
        , m_constructor(std::move(cons))
    {
        ValidateHierarchy();
    }

    const std::string& name() const noexcept
    {
        return m_name;
    }
    const MetaClass* parent() const noexcept
    {
        return mp_parent;
    }

    // ���ͼ̳й�ϵ���
    bool isKindOf(const MetaClass* pBase) const noexcept
    {
        for (const MetaClass* p = this; p != nullptr; p = p->mp_parent)
        {
            if (p == pBase)
                return true;
        }
        return false;
    }

    // ��ȫ��������
    MetaObject* create() const
    {
        if (!m_constructor)
        {
            throw std::runtime_error("Class cannot be instantiated");
        }
        return m_constructor();
    }

private:
    // ��֤��̳й�ϵ�Ϸ���
    void ValidateHierarchy() const
    {
        for (const MetaClass* p = mp_parent; p != nullptr; p = p->mp_parent)
        {
            if (p == this)
            {
                throw std::logic_error("Cyclic class hierarchy detected");
            }
        }
    }

    std::string m_name;
    const MetaClass* mp_parent;
    ObjectConstructor m_constructor;
};

//---------------------- ����Ԫ��ģ�� ----------------------
template <class T>
class ConcretMeta : public MetaClass
{
    static_assert(std::is_base_of<MetaObject, T>::value, "T must inherit from MetaObject");

public:
    static const MetaClass* instance() noexcept
    {
        static ConcretMeta s_meta;
        return &s_meta;
    }

private:
    // �Զ��Ƶ������͸�����Ϣ
    ConcretMeta()
        : MetaClass(T::StaticClassName(),                           // ����
                    ConcretMeta<typename T::MetaBase>::instance(),  // ����Ԫ��Ϣ
                    []() -> MetaObject* { return new T(); }         // ������
          )
    {
    }
};
#endif