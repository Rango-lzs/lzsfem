
#pragma once
#include "fem_export.h"
#include "FECoreFactory.h"
#include <vector>
#include <map>
#include <string.h>
#include <stdio.h>

//-----------------------------------------------------------------------------
//! This is the FECore kernel class that manages the interactions between the 
//! different modules. In particular, it manages the factory classes
//! which are responsible for the creation of different classes that are registered
//! with the kernel.
//! Used as Factory store to hold the FEFactory and 
class FEM_EXPORT FECoreKernel
{
public:
	// Do not call this function from a plugin as it will not return the correct
	// instance. Instead, use the FECoreKernel object that is passed in the PluginInitialize method
	static FECoreKernel& GetInstance();
public:
	static const char* SuperClassString(unsigned int sid);

	static std::map<unsigned int, const char*>	GetSuperClassMap();

public:
	//! Register a class with the framework
	void RegisterFactory(FECoreFactory* ptf);

	//! Create a specific using a superclass ID and an alias
	FECoreBase* Create(int superClassID, const char* szalias, FEModel* pfem);

	//! Create a class from its base class name and type string
	FECoreBase* Create(const char* baseClassName, const char* typeStr, FEModel* pfem);

	//! Create a specific class
	FECoreBase* CreateClass(const char* szclassName, FEModel* fem);

	//! Create a class from a class descriptor
	FECoreBase* Create(int superClassID, FEModel* pfem, const FEClassDescriptor& cd);

	//! Get the number of registered factory classes
	int FactoryClasses();

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int i);

	//! return a factory class
	const FECoreFactory* GetFactoryClass(int superClassID, int i);

	//! Get the index of a class factory (NOTE: this is a slow function!)
	int GetFactoryIndex(int superClassId, const char* sztype);

	//! find a factory class
	FECoreFactory* FindFactoryClass(int classID, const char* sztype);

	//! remove a factory class
	bool UnregisterFactory(FECoreFactory* ptf);

	//! unregister factories from allocator
	void UnregisterFactories(int alloc_id);

	FECoreBase* CreateInstance(const FECoreFactory* fac, FEModel* fem);

private:
	std::vector<FECoreFactory*>			m_Fac;	// list of registered factory classes
	
private: // make singleton
	FECoreKernel();
	FECoreKernel(const FECoreKernel&) = delete;
	void operator = (const FECoreKernel&) = delete;
private:
	static FECoreKernel* m_pKernel;	// the one-and-only kernel object
};

//-----------------------------------------------------------------------------
//! This class helps with the registration of a class with the framework
template <typename T> class FERegisterClass_T : public FECoreFactory
{
public:
	FERegisterClass_T(SUPER_CLASS_ID sid, const char* szclass, const char* szbase, const char* szalias, int spec = -1) : FECoreFactory(sid, szclass, szbase, szalias, spec)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterFactory(this);
	}
	FECoreBase* Create(FEModel* pfem) const { return new T(pfem); }
};

//-----------------------------------------------------------------------------
// Register a factory class
#define REGISTER_FECORE_FACTORY(theFactory) \
	static theFactory _##theFactory##_rc;

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS(theClass, ...) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theClass::superClassID(), #theClass, theClass::BaseClassName(), __VA_ARGS__);

//-----------------------------------------------------------------------------
// Register a class using default creation parameters
#define REGISTER_FECORE_CLASS_EXPLICIT(theClass, theID, ...) \
	static FERegisterClass_T<theClass> _##theClass##_rc(theID, #theClass, theClass::BaseClassName(), __VA_ARGS__);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T(theClass, theArg, theName) \
	static FERegisterClass_T<theClass<theArg> > _##theClass##theArg##_rc(theClass<theArg>::superClassID(), #theClass, theClass<theArg>::BaseClassName(), theName);

//-----------------------------------------------------------------------------
// version for classes that require template arguments
#define REGISTER_FECORE_CLASS_T2(theClass, theArg1, theArg2, theName) \
	static FERegisterClass_T<theClass<theArg1, theArg2> > _##theClass##theArg1##theArg2##_rc(theClass<theArg1, theArg2>::superClassID(), #theClass, theClass<theArg1, theArg2>::BaseClassName(), theName);

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(const char* sztype, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TBase*>(fecore.Create(TBase::superClassID(), sztype, pfem));
//	return static_cast<TBase*>(fecore.Create(TBase::BaseClassName(), sztype, pfem));
}

//-----------------------------------------------------------------------------
// Create an instance of a class.
// This assumes that TBase is derived from FECoreBase and defines a class ID. 
template <typename TBase> inline TBase* fecore_new(int classIndex, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	const FECoreFactory* f = fecore.GetFactoryClass(TBase::superClassID(), classIndex);
	if (f) return static_cast<TBase*>(f->Create(pfem));
	else return nullptr;
}

//-----------------------------------------------------------------------------
template <typename TClass> inline TClass* fecore_new_class(const char* szclass, FEModel* fem)
{
	int superId = TClass::superClassID();
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TClass*>(fecore.CreateClass(szclass, fem));
}

//-----------------------------------------------------------------------------
#define fecore_alloc(theClass, fem) fecore_new_class<theClass>(#theClass, fem)

//-----------------------------------------------------------------------------
// Three-parameter form of the fecore_new function for situations where the base class does not 
// define the classID value.
template <typename TBase> inline TBase* fecore_new(int sid, const char* sztype, FEModel* pfem)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return static_cast<TBase*>(fecore.Create(sid, sztype, pfem));
}

//=============================================================================
// TODO: Move all this stuff to sdk.h

//-----------------------------------------------------------------------------
// Template class for factory classes for plugins
template <typename T, SUPER_CLASS_ID sid> class FEPluginFactory_T : public FECoreFactory
{
public:
	FEPluginFactory_T(const char* sz) : FECoreFactory(sid, nullptr, sz, nullptr){}
	FECoreBase* Create(FEModel* pfem) const { return new T(pfem); }
};

//------------------------------------------------------------------------------
// This is for functions exported from a plugin
#ifdef WIN32
#define FECORE_EXPORT extern "C" __declspec(dllexport)
#else
#define FECORE_EXPORT extern "C"
#endif
