/*****************************************************************//**
 * \file   meta_class.cpp
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#include "femcore/RTTI/MetaClass.h"

 MetaClass::MetaClass(std::string name, const MetaClass* pParent, ObjectConstructor cons)
	 :m_name(name)
	 ,mp_parent(pParent)
	 ,m_constructor(std::move(cons))
{

}

const std::string& MetaClass::name() const
{
	return m_name;
}

const MetaClass* MetaClass::parent() const
{
	return mp_parent;
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
