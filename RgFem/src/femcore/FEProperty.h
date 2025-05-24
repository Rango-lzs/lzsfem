#pragma once
#include <vector>
#include "FEM_EXPORT.h"
#include "fecore_enum.h"

//-----------------------------------------------------------------------------
class FEObjectBase;
class DumpStream;

//-----------------------------------------------------------------------------
//! A property of a class reflects a member variable of the class that is a 
//! pointer to a FEObjectBase derived class. 
class FEM_EXPORT FEProperty
{
public:
	enum Flags
	{
		Optional		= 0x00,
		Required		= 0x01,		// the property is required (default)
		Preferred		= 0x02,		// the property is not required, but a default should be allocated when possible.
		Reference       = 0x04,		// references another class in the model
		Fixed			= 0x08,		// fixed properties are fixed type class members
		TopLevel		= 0x10,		// This is a "top-level" property. 
	};

private:
	//! Name of the property.
	//! Note that the name is not copied so it must point to a static string.
	const char*		m_szname;
	const char*		m_szlongname;	// long name (optional; used in FEBio Studio)
	unsigned int	m_flags;		// bitwise or of flags defined above 

	const char* m_szdefaultType;	// default type string (used by FEBio Studio to initialize required properties).

protected:
	const char* m_className;	// name of class that can be assigned to this

public:
	// Set\Get the name of the property
	FEProperty& SetName(const char* sz);
	const char* GetName() const;

	// get\set the long name
	FEProperty& SetLongName(const char* sz);
	const char* GetLongName() const;

	// get the class name
	const char* GetClassName() const { return m_className; }

	// is the property required
	bool IsRequired() const { return (m_flags & Required) != 0; }

	// is the property preferred
	bool IsPreferred() const { return (m_flags & Preferred) != 0; }

	// is this a reference property
	bool IsReference() const { return (m_flags & Reference) != 0; }

	// is this a top-level property
	bool IsTopLevel() const { return (m_flags & TopLevel) != 0; }

	// set the flags
	void SetFlags(unsigned int flags) { m_flags = flags; }

	// add a flag
	void AddFlag(unsigned int flag) { m_flags |= flag; }

	// get the flags
	unsigned int Flags() const { return m_flags; }

	// get default type (can be null)
	const char* GetDefaultType() const;

	// set the default type
	FEProperty& SetDefaultType(const char* szdefType);

public: // these functions have to be implemented by derived classes

	//! helper function for identifying if this is an array property or not
	virtual bool IsArray() const = 0;

	//! see if the pc parameter is of the correct type for this property
	virtual bool IsType(FEObjectBase* pc) const = 0;

	//! set the property
	virtual void SetProperty(FEObjectBase* pc) = 0;

	//! return the size of the property
	virtual int size() const = 0;

	//! return a specific property by index
	virtual FEObjectBase* get(int i) = 0;

	//! return a specific property by name
	virtual FEObjectBase* get(const char* szname) = 0;

	//! return a specific property by ID
	virtual FEObjectBase* getFromID(int nid) = 0;

	//! serialize property data
	virtual void Serialize(DumpStream& ar) = 0;

	//! initializatoin
	virtual bool Init() = 0;

	//! validation
	virtual bool Validate() = 0;

	//! Get the parent of this property
	FEObjectBase* GetParent() { return m_pParent; }

	//! Set the parent of this property
	virtual void SetParent(FEObjectBase* parent) { m_pParent = parent; }

	//! Get the class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_superClassID; }

public:
	virtual ~FEProperty();

protected:
	//! some helper functions for reading, writing properties
	void Write(DumpStream& ar, FEObjectBase* pc);
	FEObjectBase* Read(DumpStream& ar);

protected:
	// This class should not be created directly
	FEProperty(SUPER_CLASS_ID classID);

protected:
	FEObjectBase*		m_pParent;		//!< pointer to the parent class (i.e. the class that defines this property)
	SUPER_CLASS_ID	m_superClassID;	//!< The super class ID
};
