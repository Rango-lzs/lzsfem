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
#include <string>

class MetaClass;
class FEM_EXPORT MetaObject
{
public:
	using BaseClass = void;
	virtual ~MetaObject() = 0;

	virtual const MetaClass* meta() const = 0;  //called by instance
	static const MetaClass* staic_meta();  //called by class
	static MetaObject* meta_cast(MetaObject* pOther);
    static std::string class_name();
	bool isKindOf(const MetaClass* pMeta) const;

public:
	MetaObject() = default;
};

#endif