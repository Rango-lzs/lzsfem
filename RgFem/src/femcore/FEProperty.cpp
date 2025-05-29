#include "FEProperty.h"
#include "basicio/DumpStream.h"
#include "FEObjectBase.h"
//-----------------------------------------------------------------------------
FEProperty::FEProperty() : m_szname(nullptr), m_className(nullptr), m_szlongname(nullptr), m_szdefaultType(nullptr), m_flags(0), m_superClassID() {}

//-----------------------------------------------------------------------------
FEProperty::~FEProperty(){}

//-----------------------------------------------------------------------------
//! Set the name of the property.
//! Note that the name is not copied so it must point to a static string.
FEProperty& FEProperty::SetName(const char* sz)
{
	m_szname = sz;
	return *this;
}

//-----------------------------------------------------------------------------
//! Return the name of this property
const char* FEProperty::GetName() const { return m_szname; }

//-----------------------------------------------------------------------------
FEProperty& FEProperty::SetLongName(const char* sz)
{
	m_szlongname = sz;
	return *this;
}

//-----------------------------------------------------------------------------
const char* FEProperty::GetLongName() const { return m_szlongname; }

//-----------------------------------------------------------------------------
const char* FEProperty::GetDefaultType() const
{
	return m_szdefaultType;
}

//-----------------------------------------------------------------------------
FEProperty& FEProperty::SetDefaultType(const char* szdefType)
{
	m_szdefaultType = szdefType;
	return *this;
}

//-----------------------------------------------------------------------------
void FEProperty::Write(DumpStream& ar, FEObjectBase* pc)
{
	int nflag = (pc == 0 ? 0 : 1);
	ar << nflag;
	if (nflag)
	{		
		ar << pc->GetTypeStr();	
		pc->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
FEObjectBase* FEProperty::Read(DumpStream& ar)
{
	int nflag = 0;
	FEObjectBase* pm = 0;
	ar >> nflag;
	if (nflag)
	{
		char sz[256];
		ar >> sz;
		pm = RANGO_NEW<FEObjectBase>(&ar.GetFEModel(), sz);	
		pm->Serialize(ar);

		// TODO: Do I really need to do this here?
		//pm->Init();
	}
	return pm;
}
