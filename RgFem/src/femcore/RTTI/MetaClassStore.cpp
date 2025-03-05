/*********************************************************************
 * \file   MetaClassStore.h
 * \brief
 *
 * \author Leizs
 * \date   March 2025
 *********************************************************************/
#include "MetaClassStore.h"
#include "MetaClass.h"

const MetaClass* MetaClassStore::get(const std::string& rxName) const
{
    return m_meta_store.at(rxName);
}
void MetaClassStore::insert(const MetaClass* pMeta)
{
    m_meta_store.emplace(pMeta->name(), pMeta);
}

MetaClassStore::MetaClassStore()
{

}
