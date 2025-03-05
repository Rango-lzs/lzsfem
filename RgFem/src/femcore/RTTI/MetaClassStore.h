/*********************************************************************
 * \file   MetaClassStore.h
 * \brief  
 * 
 * \author Leizs
 * \date   March 2025
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"
#include "femcore/RTTI/MetaClass.h"
#include <map>
#include <string>

class MetaClass;

//��������ڴ洢�������MetaClass��Ϣ��
//��Ҫͨ������������ȡ���Ӧ��MetaClass���Ӷ��������ʵ��
//�����Ӧ���Ǹ�����
class FEM_EXPORT MetaClassStore
{
public:
    static MetaClassStore* instance()
    {
        static MetaClassStore metaStore;
        return &metaStore;
    }
    const MetaClass* get(const std::string& rxName) const
    {
        return m_meta_store.at(rxName);
    }
    void insert(const MetaClass* pMeta)
    {
        m_meta_store.emplace(pMeta->name(), pMeta);
    }

private:
    MetaClassStore();
    std::map<std::string, const MetaClass*> m_meta_store;
};
