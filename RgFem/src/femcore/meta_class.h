/*****************************************************************//**
 * \file   meta_class.h
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/
#pragma once
#ifndef META_CLASS_H
#define META_CLASS_H

using ObjectConstructor = std::function<MetaObject* ()>;

class MetaClass
{
public:
	//��Ҫͨ�����캯����������MetaClass���ܹ�ʵ����MetaObject�Ĺ��캯��
	MetaClass(const std::string& name, const MetaClass* pParent, ObjectConstructor cons);
	~MetaClass() = default;

	std::string name();
	const MetaClass* parent();
	bool isKindOf(const MetaClass* pParent);

	MetaObject* create();

private:
	std::string m_name;
	const MetaClass* mp_parent;
	ObjectConstructor m_constructor;
};