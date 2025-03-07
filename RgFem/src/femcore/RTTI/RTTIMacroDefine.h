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
    static const MetaClass* staic_meta();                                                                              \
    static std::string class_name();                                                                                   \
    static MetaObject* meta_cast(MetaObject* pOther)


#define DEFINE_META_CLASS(DERIVE_CLASS, BASE_CLASS)                                                                    \
    const const MetaClass* DERIVE_CLASS::meta() const                                                                  \
    {                                                                                                                  \
        return DERIVE_CLASS::staic_meta();                                                                             \
    }                                                                                                                  \
    const MetaClass* DERIVE_CLASS::staic_meta()                                                                        \
    {                                                                                                                  \
        return ConcretMeta<DERIVE_CLASS>::instance();                                                                  \
    }                                                                                                                  \
    std::string DERIVE_CLASS::class_name()                                                                             \
    {                                                                                                                  \
        return #DERIVE_CLASS;                                                                                          \
    }                                                                                                                  \
    MetaObject* DERIVE_CLASS::meta_cast(MetaObject* pOther)                                                            \
    {                                                                                                                  \
        if (!pOther)                                                                                                   \
            return nullptr;                                                                                            \
        return pOther->isKindOf(MetaObject::staic_meta()) ? static_cast<MetaObject*>(pOther) : nullptr;                \
    }                                                                                                                  \
    static const MetaClass* s_pMeta = ConcretMeta<DERIVE_CLASS>::instance()