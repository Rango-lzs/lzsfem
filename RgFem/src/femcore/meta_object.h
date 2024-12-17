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

	virtual const MetaClass* meta();
	static const MetaClass* staic_meta();
	bool isKindOf(const MetaClass* pMeta);

protected:
	MetaObject() = default;
};

#endif