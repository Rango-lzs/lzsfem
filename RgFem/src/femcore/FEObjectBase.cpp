
#include "FEObjectBase.h"
#include "FEModelParam.h"
#include "logger/log.h"
#include <sstream>

//-----------------------------------------------------------------------------
FEObjectBase::FEObjectBase(FEModel* fem) : m_fem(fem)
{ 
	m_nID = -1;
}

//-----------------------------------------------------------------------------
//! destructor does nothing for now.
FEObjectBase::~FEObjectBase()
{
}


//-----------------------------------------------------------------------------
//! Sets the user defined name of the component
void FEObjectBase::SetName(const std::string& name)
{ 
	m_name = name;
}

//-----------------------------------------------------------------------------
//! Return the name
const std::string& FEObjectBase::GetName() const
{ 
	return m_name; 
}


//-----------------------------------------------------------------------------
//! return the component ID
int FEObjectBase::GetID() const { return m_nID; }

//-----------------------------------------------------------------------------
//! set the component ID
void FEObjectBase::SetID(int nid) { m_nID = nid; }

//-----------------------------------------------------------------------------
//! Get the FE model
FEModel* FEObjectBase::GetFEModel() const { return m_fem; }

//-----------------------------------------------------------------------------
void FEObjectBase::SetFEModel(FEModel* fem) { m_fem = fem; }


//-----------------------------------------------------------------------------
void FEObjectBase::Serialize(DumpStream& ar)
{
	//// do base class first
	//FEParamContainer::Serialize(ar);

	//// serialize name
	//if (ar.IsShallow() == false)
	//{
	//	ar & m_name;
	//	ar & m_nID;
	//}

	//if (ar.IsShallow() == false) ar & m_pParent;

	//// serialize all the properties
	//int NP = (int)m_Prop.size();
	//for (int i = 0; i<NP; ++i)
	//{
	//	FEProperty* prop = m_Prop[i];
	//	prop->SetParent(this);
	//	prop->Serialize(ar);
	//}
}

//-----------------------------------------------------------------------------
void FEObjectBase::SaveClass(DumpStream& ar, FEObjectBase* a)
{
	/*assert(ar.IsSaving());
	int classID = 0;
	if (a == nullptr) { ar << classID; return; }
	classID = a->GetSuperClassID();
	assert(classID != FEINVALID_ID);
	const char* sztype = a->GetTypeStr();
	ar << classID;
	ar << sztype;*/
}

FEObjectBase* FEObjectBase::LoadClass(DumpStream& ar, FEObjectBase* a)
{
	//assert(ar.IsLoading());

	//int classID = 0;
	//ar >> classID;
	//if (classID == FEINVALID_ID) return nullptr;

	//char sztype[256] = { 0 };
	//ar >> sztype;

	//// instantiate the class
	//a = fecore_new<FEObjectBase>(classID, sztype, &ar.GetFEModel());
	//assert(a);

	//if (a == nullptr) throw DumpStream::ReadError();

	//return a;
}

//-----------------------------------------------------------------------------
bool FEObjectBase::Validate()
{
	// validate parameters
	//FEParameterList& pl = GetParameterList();
	//int N = pl.Parameters();
	//list<FEParam>::iterator pi = pl.first();
	//for (int i = 0; i < N; ++i, pi++)
	//{
	//	FEParam& p = *pi;
	//	if (p.is_valid() == false)
	//	{
	//		stringstream ss;
	//		ss << GetName() << "." << p.name();
	//		string paramName = ss.str();
	//		feLogError("Invalid value for parameter: %s", paramName.c_str());
	//		return false;
	//	}
	//}

	//// check properties
	//const int nprop = (int)m_Prop.size();
	//for (int i = 0; i<nprop; ++i)
	//{
	//	FEProperty* pi = m_Prop[i];
	//	if (pi)
	//	{
	//		if (pi->Validate() == false) return false;
	//	}
	//}

	return true;
}

//-----------------------------------------------------------------------------
bool FEObjectBase::Init()
{
	// call init on model parameters
	FEParameterList& PL = GetParameterList();
	FEParamIterator it = PL.first();
	for (int i = 0; i < PL.Parameters(); ++i, ++it)
	{
		FEParam& pi = *it;
		if (pi.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamDouble& pd = pi.value<FEParamDouble>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
		else if (pi.type() == FE_PARAM_VEC3D_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamVec3& pd = pi.value<FEParamVec3>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
		else if (pi.type() == FE_PARAM_MAT3D_MAPPED)
		{
			for (int j = 0; j < pi.dim(); ++j)
			{
				FEParamMat3d& pd = pi.value<FEParamMat3d>(j);
				if (pd.Init() == false)
				{
					feLogError("Failed to initialize parameter %s", pi.name());
					return false;
				}
			}
		}
	}

	// check the parameter ranges
	if (Validate() == false) return false;

	// initialize properties
	/*const int nprop = (int)m_Prop.size();
	for (int i = 0; i<nprop; ++i)
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Init() == false)
			{
				feLogError("The property \"%s\" failed to initialize or was not defined.", pi->GetName());
				return false;
			}
		}
		else
		{
			feLogError("A nullptr was set for property i");
			return false;
		}
	}*/
	return true;
}

//-----------------------------------------------------------------------------
void FEObjectBase::AddProperty(FEProperty* pp, const char* sz, unsigned int flags)
{
	pp->SetName(sz);
	pp->SetLongName(sz);
	pp->SetFlags(flags);
	pp->SetParent(this);
	m_Prop.push_back(pp);
}

//-----------------------------------------------------------------------------
void FEObjectBase::RemoveProperty(int i)
{
	m_Prop[i] = nullptr;
}

//-----------------------------------------------------------------------------
void FEObjectBase::ClearProperties()
{
	for (int i = 0; i < m_Prop.size(); ++i)
	{
		delete m_Prop[i];
	}
	m_Prop.clear();
}

//-----------------------------------------------------------------------------
int FEObjectBase::Properties()
{
	int N = (int)m_Prop.size();
	int n = 0;
	for (int i = 0; i<N; ++i) n += m_Prop[i]->size();
	return n;
}

//-----------------------------------------------------------------------------
int FEObjectBase::FindPropertyIndex(const char* sz)
{
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		const FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
FEProperty* FEObjectBase::FindProperty(const char* sz, bool searchChildren)
{
	// first, search the class' properties
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return pm;
	}

	// the property, wasn't found so look into the properties' properties
	if (searchChildren)
	{
		for (int i = 0; i < NP; ++i)
		{
			FEProperty* pm = m_Prop[i];
			if (pm)
			{
				int m = pm->size();
				for (int j = 0; j < m; ++j)
				{
					FEObjectBase* pcj = pm->get(j);
					if (pcj)
					{
						// Note: we don't search children's children!
						FEProperty* pj = pcj->FindProperty(sz);
						if (pj) return pj;
					}
				}
			}
		}
	}

	return nullptr;
}

//-----------------------------------------------------------------------------
FEObjectBase* FEObjectBase::GetProperty(int n)
{
	int N = (int)m_Prop.size();
	int m = 0;
	for (int i = 0; i<N; ++i)
	{
		FEProperty* pm = m_Prop[i];
		int l = pm->size();
		if (m + l > n) return pm->get(n - m);
		m += l;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool FEObjectBase::SetProperty(int i, FEObjectBase* pb)
{
	FEProperty* pm = m_Prop[i];
	if (pm->IsType(pb))
	{
		pm->SetProperty(pb);
		if (pb) pb->SetParent(this);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Set a property via name
bool FEObjectBase::SetProperty(const char* sz, FEObjectBase* pb)
{
	FEProperty* prop = FindProperty(sz);
	if (prop == nullptr) return false;

	if (prop->IsType(pb))
	{
		prop->SetProperty(pb);
		if (pb) pb->SetParent(this);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! number of parameters
int FEObjectBase::Parameters() const
{
	return GetParameterList().Parameters();
}

//-----------------------------------------------------------------------------
FEParam* FEObjectBase::FindParameter(const ParamString& s)
{
	// first search the parameter list
	FEParam* p = FEParamContainer::FindParameter(s);
	if (p) return p;

	// next, let's try the property list
	int NP = (int)m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		// get the property
		FEProperty* mp = m_Prop[i];

		// see if matches
		if (s == mp->GetName())
		{
			if (mp->IsArray())
			{
				// get the number of items in this property
				int nsize = mp->size();
				int index = s.Index();
				if ((index >= 0) && (index < nsize))
				{
					return mp->get(index)->FindParameter(s.next());
				}
				else
				{
					int nid = s.ID();
					if (nid != -1)
					{
						FEObjectBase* pc = mp->getFromID(nid);
						if (pc) return pc->FindParameter(s.next());
					}
					else if (s.IDString())
					{
						FEObjectBase* c = mp->get(s.IDString());
						if (c) return c->FindParameter(s.next());
					}
				}
			}
			else
			{
				FEObjectBase* pc = mp->get(0);
				return (pc ? pc->FindParameter(s.next()) : nullptr);
			}
		}
	}

	return nullptr;
}

//-----------------------------------------------------------------------------
//! return the property (or this) that owns a parameter
FEObjectBase* FEObjectBase::FindParameterOwner(void* pd)
{
	// see if this class is the owner of the data pointer
	FEParam* p = FindParameterFromData(pd);
	if (p) return this;

	// it's not se let's check the properties
	int NP = PropertyClasses();
	for (int i = 0; i < NP; ++i)
	{
		FEProperty* pi = PropertyClass(i);
		int n = pi->size();
		for (int j = 0; j < n; ++j)
		{
			FEObjectBase* pcj = pi->get(j);
			if (pcj)
			{
				FEObjectBase* pc = pcj->FindParameterOwner(pd);
				if (pc) return pc;
			}
		}
	}

	// sorry, no luck
	return nullptr;
}

//-----------------------------------------------------------------------------
//! return the number of properties defined
int FEObjectBase::PropertyClasses() const 
{ 
	return (int)m_Prop.size(); 
}

//! return a property
FEProperty* FEObjectBase::PropertyClass(int i)
{ 
	return m_Prop[i]; 
}

//-----------------------------------------------------------------------------
FEProperty* FEObjectBase::FindProperty(const ParamString& prop)
{
	int NP = (int)m_Prop.size();
	for (int i = 0; i < NP; ++i)
	{
		FEProperty* mp = m_Prop[i];

		if (prop == mp->GetName())
		{
			if (mp->IsArray())
			{
				// get the number of items in this property
				int nsize = mp->size();
				int index = prop.Index();
				if ((index >= 0) && (index < nsize))
				{
					FEObjectBase* pc = mp->get(index);
					if (pc)
					{
						ParamString next = prop.next();
						if (next.count() == 0) return nullptr;
						else return pc->FindProperty(next);
					}
				}
				else
				{
					int nid = prop.ID();
					if (nid != -1)
					{
						FEObjectBase* pc = mp->getFromID(nid);
						// TODO: What to do here? 
						assert(false);
					}
					else if (prop.IDString())
					{
						FEObjectBase* pc = mp->get(prop.IDString());
						if (pc)
						{
							ParamString next = prop.next();
							if (next.count() == 0) return nullptr;
							else return pc->FindProperty(next);
						}
					}
				}
			}
			else
			{
				ParamString next = prop.next();
				if (next.count() == 0) return mp;
				else return FindProperty(next);
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
FEObjectBase* FEObjectBase::GetProperty(const ParamString& prop)
{
	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{	
		FEProperty* mp = m_Prop[i];

		if (prop == mp->GetName())
		{
			if (mp->IsArray())
			{
				// get the number of items in this property
				int nsize = mp->size();
				int index = prop.Index();
				if ((index >= 0) && (index < nsize))
				{
					FEObjectBase* pc = mp->get(index);
					if (pc)
					{
						ParamString next = prop.next();
						if (next.count() == 0) return pc;
						else return pc->GetProperty(next);
					}
				}
				else
				{
					int nid = prop.ID();
					if (nid != -1)
					{
						FEObjectBase* pc = mp->getFromID(nid);
					}
					else if (prop.IDString())
					{
						FEObjectBase* pc = mp->get(prop.IDString());
						if (pc)
						{
							ParamString next = prop.next();
							if (next.count() == 0) return pc;
							else return pc->GetProperty(next);
						}
					}
				}
			}
			else
			{
				FEObjectBase* pc = mp->get(0);
				ParamString next = prop.next();
				if (next.count() == 0) return pc;
				else return pc->GetProperty(next);
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
bool FEObjectBase::BuildClass()
{
	GetParameterList();

	for (int i = 0; i < PropertyClasses(); ++i)
	{
		FEProperty* pp = PropertyClass(i);
		int m = pp->size();
		for (int j = 0; j < m; ++j)
		{
			FEObjectBase* pj = pp->get(j);
			if (pj) pj->BuildClass();
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEObjectBase::UpdateParams()
{
	return true;
}
