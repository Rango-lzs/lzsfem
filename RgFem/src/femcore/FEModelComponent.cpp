#include "FEModelComponent.h"
#include "femcore/FEModel.h"
#include "basicio/DumpStream.h"
#include <string.h>

DEFINE_META_CLASS(FEModelComponent, FEObjectBase, "");

//-----------------------------------------------------------------------------
FEModelComponent::FEModelComponent() : FEObjectBase()
{
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{
	
}

//-----------------------------------------------------------------------------
void FEModelComponent::Update()
{

}

//-----------------------------------------------------------------------------
double FEModelComponent::CurrentTime() const
{
	return GetFEModel()->GetTime().currentTime;
}

//-----------------------------------------------------------------------------
double FEModelComponent::CurrentTimeIncrement() const
{
	return GetFEModel()->GetTime().timeIncrement;
}

//-----------------------------------------------------------------------------
double FEModelComponent::GetGlobalConstant(const char* sz) const
{
	return GetFEModel()->GetGlobalConstant(sz);
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetDOFIndex(const char* szvar, int n) const
{
	return GetFEModel()->GetDOFIndex(szvar, n);
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetDOFIndex(const char* szdof) const
{
	return GetFEModel()->GetDOFIndex(szdof);
}

//-----------------------------------------------------------------------------
//! Get the model's mesh
FEMesh& FEModelComponent::GetMesh()
{
	return GetFEModel()->GetMesh();
}

//-----------------------------------------------------------------------------
const FETimeInfo& FEModelComponent::GetTimeInfo() const
{
	return GetFEModel()->GetTime();
}

//-----------------------------------------------------------------------------
void FEModelComponent::AttachLoadController(const char* szparam, int lc)
{
	FEParam* p = GetParameter(szparam); assert(p);
	if (p)
	{
		FEModel* fem = GetFEModel();
		fem->AttachLoadController(p, lc);
	}
}

//-----------------------------------------------------------------------------
void FEModelComponent::AttachLoadController(void* pd, int lc)
{
	FEParam* p = FindParameterFromData(pd); assert(p);
	if (p)
	{
		FEModel* fem = GetFEModel();
		fem->AttachLoadController(p, lc);
	}
}
