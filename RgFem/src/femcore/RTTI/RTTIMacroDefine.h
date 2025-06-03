/*********************************************************************
 * \file   RTTIMacroDefine.h
 * \brief
 *
 * \author Leizs
 * \date   March 2025
 *********************************************************************/

#pragma once

#define DECLARE_META_CLASS(DERIVE_CLASS, BASE_CLASS)                                                                   \
public:                                                                                                                \
    using BaseClass = BASE_CLASS;                                                                                      \
    virtual const const MetaClass* meta() const;                                                                       \
    static const MetaClass* static_meta();                                                                             \
    static std::string class_name();                                                                                   \
    static std::string alias_name();                                                                                   \
    static MetaObject* meta_cast(MetaObject* pOther)


#define DEFINE_META_CLASS(DERIVE_CLASS, BASE_CLASS, ALIAS_NAME)                                                        \
    const const MetaClass* DERIVE_CLASS::meta() const                                                                  \
    {                                                                                                                  \
        return DERIVE_CLASS::static_meta();                                                                            \
    }                                                                                                                  \
    const MetaClass* DERIVE_CLASS::static_meta()                                                                       \
    {                                                                                                                  \
        return ConcreteMeta<DERIVE_CLASS>::instance();                                                                 \
    }                                                                                                                  \
    std::string DERIVE_CLASS::class_name()                                                                             \
    {                                                                                                                  \
        return #DERIVE_CLASS;                                                                                          \
    }                                                                                                                  \
    std::string DERIVE_CLASS::alias_name()                                                                             \
    {                                                                                                                  \
        return ALIAS_NAME;                                                                                             \
    }                                                                                                                  \
    MetaObject* DERIVE_CLASS::meta_cast(MetaObject* pOther)                                                            \
    {                                                                                                                  \
        if (!pOther)                                                                                                   \
            return nullptr;                                                                                            \
        return pOther->isKindOf(MetaObject::static_meta()) ? static_cast<MetaObject*>(pOther) : nullptr;               \
    }                                                                                                                  \
    static const MetaClass* DERIVE_CLASS##_pMeta = ConcreteMeta<DERIVE_CLASS>::instance()