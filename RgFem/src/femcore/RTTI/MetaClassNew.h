#pragma once

#ifndef META_CLASS_H
#define META_CLASS_H

#include <functional>
#include <string>
#include <type_traits>

class MetaObject;

//---------------------- 类型别名声明 ----------------------
using ObjectConstructor = std::function<MetaObject*()>;

//---------------------- 宏定义（简化用户声明）----------------------
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

//---------------------- MetaClass 核心定义 ----------------------
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

    // 类型继承关系检查
    bool isKindOf(const MetaClass* pBase) const noexcept
    {
        for (const MetaClass* p = this; p != nullptr; p = p->mp_parent)
        {
            if (p == pBase)
                return true;
        }
        return false;
    }

    // 安全创建对象
    MetaObject* create() const
    {
        if (!m_constructor)
        {
            throw std::runtime_error("Class cannot be instantiated");
        }
        return m_constructor();
    }

private:
    // 验证类继承关系合法性
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

//---------------------- 具体元类模板 ----------------------
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
    // 自动推导类名和父类信息
    ConcretMeta()
        : MetaClass(T::StaticClassName(),                           // 类名
                    ConcretMeta<typename T::MetaBase>::instance(),  // 父类元信息
                    []() -> MetaObject* { return new T(); }         // 构造器
          )
    {
    }
};
#endif