/*****************************************************************//**
 * \file   meta_object.h
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once
#ifndef META_OBJECT_H
#define META_OBJECT_H

#include "femcore/fem_export.h"

class MetaClass;
class FEM_EXPORT MetaObject
{
public:
	virtual ~MetaObject() = default;

	virtual const MetaClass* meta() const;  //called by instance
	static const MetaClass* staic_meta();  //called by class
	static MetaObject* meta_cast(MetaObject* pOther);	
	bool isKindOf(const MetaClass* pMeta) const;

protected:
	MetaObject() = default;
};

class DemoObject : public MetaObject
{
public:
	DemoObject() {}

	virtual const const MetaClass* meta()
	{
		return mpMeta;
	}

	static const MetaClass* staic_meta()
	{
		return mpMeta;
	}
public:
	void member()
	{

	}
private:
	static MetaClass* mpMeta;
};

//MetaClass* DemoObject::mpMeta = new MetaClass("DemoObject", MetaObject::meta(), []() { return new DemoObject(); });

#endif