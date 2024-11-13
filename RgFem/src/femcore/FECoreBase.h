/*****************************************************************//**
 * \file   FECoreObject.h
 * \brief  
 * 
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#pragma once
#include "FEParameterList.h"
#include "fecore_enum.h"
#include "FEProperty.h"
#include "ClassDescriptor.h"
#include <string>

//-----------------------------------------------------------------------------
class FECoreFactory;
class FEModel;
class FEStream;
class FEMetalClass;

/*
* Propertyr 和 Paramater的区别
*/

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all 
//! classes that can be registered with the framework.
class FEM_EXPORT FECoreObject
{
public:
	//! constructor
	FECoreObject(FEModel* fem);

	//! destructor
	virtual ~FECoreObject();

	//! data serialization
	void Serialize(FEStream& ar) override;

	//! Initialization
	virtual bool Init();

	//! validates all properties and parameters
	bool Validate() override;

	//! call this after the parameters are changed
	virtual bool UpdateParams();

	// Build the class' parameter and property list
	bool BuildClass();

	//! number of parameters
	int Parameters() const;

	//! return a parameter
	virtual FEParam* FindParameter(const ParamString& s) override;

	//! return the property (or this) that owns a parameter
	FECoreObject* FindParameterOwner(void* pd);

public: // interface for getting/setting properties

	//! get the number of properties
	int Properties();

	//! Set a property
	bool SetProperty(int nid, FECoreObject* pm);

	//! Set a property via name
	bool SetProperty(const char* sz, FECoreObject* pm);

	//! return a property
	virtual FECoreObject* GetProperty(int i);

	//! find a property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! return a property (class)
	FEProperty* FindProperty(const char* sz, bool searchChildren = false);

	//! return a property from a paramstring
	FECoreObject* GetProperty(const ParamString& prop);

	//! return the number of properties defined
	int PropertyClasses() const;

	//! return a property
	FEProperty* PropertyClass(int i);

	//! return the component ID
	int GetID() const;

	//! set the component ID
	void SetID(int nid);

	//! Get the FE model
	FEModel* GetFEModel() const;

	//! set the FEModel of this class (use with caution!)
	void SetFEModel(FEModel* fem);

	static void SaveClass(FEStream& ar, FECoreObject* p);
	static FECoreObject* LoadClass(FEStream& ar, FECoreObject* p);

	// set parameters through a class descriptor
	bool SetParameters(const FEClassDescriptor& cd);
	bool SetParameters(const FEClassDescriptor::ClassVariable& cv);

public:
	//! Add a property
	//! Call this in the constructor of derived classes to 
	//! build the property list
	void AddProperty(FEProperty* pp, const char* sz, unsigned int flags = FEProperty::Required);

	void RemoveProperty(int i);

	void ClearProperties();

public:
	template <class T> T* ExtractProperty(bool extractSelf = true);

private:
	std::string		m_name;			//!< user defined name of component
	FECoreObject*		m_pParent;		//!< pointer to "parent" object (if any) (NOTE: only used by materials)
	FEModel*		m_fem;			//!< the model this class belongs to

	vector<FEProperty*>	m_Prop;		//!< list of properties

private:
	int		m_nID;			//!< component ID
	const FECoreFactory*	m_fac;	//!< factory class that instantiated this class
	friend class FECoreFactory;
};

// include template property definitions
#include "FEPropertyT.h"

template <class T>	FEProperty* AddClassProperty(FECoreObject* pc, T* pp, const char* sz)
{
	FEFixedPropertyT<T>* prop = new FEFixedPropertyT<T>(pp);
	prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, FEProperty::Fixed);
	return prop;
}

template <class T> FEProperty* AddClassProperty(FECoreObject* pc, T** pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEPropertyT<T>* prop = new FEPropertyT<T>(pp);
	if (prop->GetSuperClassID() == FECLASS_ID) prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, flags);
	return prop;
}

template <class T>	FEProperty* AddClassProperty(FECoreObject* pc, std::vector<T*>* pp, const char* sz, unsigned int flags = FEProperty::Required)
{
	FEVecPropertyT<T>* prop = new FEVecPropertyT<T>(pp);
	if (prop->GetSuperClassID() == FECLASS_ID) prop->SetDefaultType(sz);
	pc->AddProperty(prop, sz, flags);
	return prop;
}

#define ADD_PROPERTY(theProp, ...) AddClassProperty(this, &theProp, __VA_ARGS__)

#define FECORE_SUPER_CLASS(a) public: static SUPER_CLASS_ID superClassID() { return a; }
//#define REGISTER_SUPER_CLASS(theClass, a) SUPER_CLASS_ID theClass::superClassID() { return a;}

template <class T> T* FECoreObject::ExtractProperty(bool extractSelf)
{
	if (extractSelf)
	{
		if (dynamic_cast<T*>(this)) return dynamic_cast<T*>(this);
	}

	int NC = Properties();
	for (int i = 0; i < NC; i++)
	{
		FECoreObject* pci = GetProperty(i);
		if (pci)
		{
			T* pc = pci->ExtractProperty<T>();
			if (pc) return pc;
		}
	}
	return nullptr;
}
