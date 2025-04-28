/*****************************************************************//**
 * \file   meta_class.cpp
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#include "femcore/RTTI/MetaClass.h"

 MetaClass::MetaClass(const std::string& className, const std::string& aliasName, MetaClass* pParent, ObjectConstructor cons)
	 :m_class_name(className)
	 ,m_alias_name(aliasName)
	 ,mp_parent(pParent)
	 ,m_constructor(std::move(cons))
{
	 pParent->m_childs.push_back(this);
}

const std::string& MetaClass::name() const
{
	return m_class_name;
}

const std::string& MetaClass::alias_name() const
{
    return m_alias_name;
}

const MetaClass* MetaClass::parent() const
{
	return mp_parent;
}

const std::list<MetaClass*>& MetaClass::childs() const
{
    return m_childs;
}

bool MetaClass::isKindOf(const MetaClass* pParent) const
{
	const MetaClass* pCurMeta = this;
	while (pCurMeta)
	{
		if (pCurMeta == pParent)
		{
			return true;
		}
		pCurMeta = pCurMeta->parent();
	}
	return false;
}

MetaObject* MetaClass::create() const
{
	return m_constructor();
}
