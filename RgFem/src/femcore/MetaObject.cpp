/*****************************************************************
 * \file   MetaObject.cpp
 * \brief
 *
 * \author 11914
 * \date   March 2025
 *********************************************************************/

#include "femcore/MetaObject.h"
#include "femcore/RTTI/MetaClass.h"

static MetaClass* gsMetaObjClass = nullptr;   //MetaObject ×ÔÉíµÄMetaClass

const MetaClass* MetaObject::meta() const
{
    return gsMetaObjClass;
}

const MetaClass* MetaObject::staic_meta()
{
    return gsMetaObjClass;
}

std::string class_name()
{
    return "MetaObject";
}

MetaObject* MetaObject::meta_cast(MetaObject* pOther)
{
    if (!pOther)
        return nullptr;
    return pOther->isKindOf(MetaObject::staic_meta()) ? static_cast<MetaObject*>(pOther) : nullptr;
}


bool MetaObject::isKindOf(const MetaClass* pMeta) const
{
    const MetaClass* pMetaThis = this->meta();
    if (pMetaThis)
    {
        return pMetaThis->isKindOf(pMeta);
    }
    return false;
}
