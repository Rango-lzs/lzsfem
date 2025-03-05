/*********************************************************************
 * \file   MetaClassStore.h
 * \brief  
 * 
 * \author Leizs
 * \date   March 2025
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"
#include <map>
#include <string>

class MetaClass;
//这个类用于存储所有类的MetaClass信息，
//需要通过类名称来获取其对应的MetaClass，从而构造类的实例
//这个类应该是个单例
class FEM_EXPORT MetaClassStore
{
public:
    static MetaClassStore* instance()
    {
        static MetaClassStore metaStore;
        return &metaStore;
    }
    const MetaClass* get(const std::string& rxName) const;
    
    void insert(const MetaClass* pMeta);
    
private:
    MetaClassStore();
    std::map<std::string, const MetaClass*> m_meta_store;
};
