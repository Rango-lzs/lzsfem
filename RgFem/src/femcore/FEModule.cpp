#include "FEModule.h"
#include <vector>
#include <cstring>

class FEModule::Impl
{
public:
	const char*			szname = 0;		// name of module
	const char*			szdesc = 0;		// description of module (optional, can be null)
	unsigned int		id = 0;			// unqiue ID (starting at one)
	int					alloc_id = -1;	// ID of allocator
	int					m_status = FEModule::RELEASED;	// Status of module
	std::vector<int>	depMods;	// module dependencies

public:
	void AddDependency(int mid)
	{
		for (size_t i = 0; i < depMods.size(); ++i)
		{
			if (depMods[i] == mid) return;
		}
		depMods.push_back(mid);
	}

	void AddDependencies(const std::vector<int>& mid)
	{
		for (size_t i = 0; i < mid.size(); ++i)
		{
			AddDependency(mid[i]);
		}
	}
};

FEModule::FEModule() : im(new FEModule::Impl)
{

}

FEModule::FEModule(const char* szname, const char* szdescription) : im(new FEModule::Impl)
{
	SetName(szname);
	SetDescription(szdescription);
}

FEModule::~FEModule()
{
	delete im;
}

int FEModule::GetModuleID() const
{
	return im->id;
}

const char* FEModule::GetName() const
{
	return im->szname;
}

const char* FEModule::GetDescription() const
{
	return im->szdesc;
}

void FEModule::SetStatus(FEModule::Status status)
{
	im->m_status = status;
}

int FEModule::GetStatus() const
{
	return im->m_status;
}

bool FEModule::HasDependent(int modId) const
{
	if (modId == im->id) return true;
	for (int i : im->depMods)
	{
		if (i == modId) return true;
	}
	return false;
}

void FEModule::AddDependency(FEModule& mod)
{
	im->AddDependency(mod.GetModuleID());
	im->AddDependencies(mod.im->depMods);
}

void FEModule::ClearDependencies()
{
	im->depMods.clear();
}

std::vector<int> FEModule::GetDependencies() const
{
	return im->depMods;
}

void FEModule::SetID(int newId)
{
	im->id = newId;
}

void FEModule::SetName(const char* szname)
{
	im->szname = szname;
}

void FEModule::SetDescription(const char* szdesc)
{
	im->szdesc = szdesc;
}
