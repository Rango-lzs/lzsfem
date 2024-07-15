
#pragma once

//-----------------------------------------------------------------------------
//! Forward declaration of the FEModel class. All classes that register
//! with the framework take a pointer to FEModel as their constructor parameter.
class FEModel;
class FECoreBase;
enum class SUPER_CLASS_ID;

//-----------------------------------------------------------------------------
//! The factory class contains the mechanism for instantiating a class.
//! Also play as the role of a RxClass 
class FEM_EXPORT FECoreFactory
{
public:
	//! constructor
	FECoreFactory(SUPER_CLASS_ID scid, const char* szclass, const char* szbase, const char* szalias, int nspec = -1);

	//! virtual constructor
	virtual ~FECoreFactory();

	//! This is the function that the kernel will use to intantiate an object
	FECoreBase* CreateInstance(FEModel* pfem) const;

public:
	// return the class name
	const char* GetClassName() const { return m_szclass; }

	// return the base class name
	const char* GetBaseClassName() const { return m_szbase; }

	// return the type string identifier
	const char* GetTypeStr() const { return m_szalias; }

	//! return the super-class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_scid; }

	//! return the module name
	unsigned int GetModuleID() const { return m_module; }

	//! set the module name
	void SetModuleID(unsigned int nid);

	//! Get the spec number
	int GetSpecID() const { return m_spec; }

	//! Set the allocator ID
	void SetAllocatorID(int alloc) { m_alloc_id = alloc; }

	//! Get the allocator ID
	int GetAllocatorID() const { return m_alloc_id; }
	
public:
	//! derived classes implement this to create an instance of a class
	virtual FECoreBase* Create(FEModel*) const = 0;

private:
	const char*		m_szclass;	//!< class name
	const char*		m_szbase;	//!< base class name
	const char*		m_szalias;	//!< class alias string
	int				m_spec;		//!< The max spec number for which this feature is defined (-1 is don't care)
	unsigned int	m_module;	//!< ID of module this class belongs to
	SUPER_CLASS_ID	m_scid;		//!< the super-class ID
	int				m_alloc_id;	//!< allocator ID
};
