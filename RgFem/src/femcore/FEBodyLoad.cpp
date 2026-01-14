
#include "FEBodyLoad.h"
#include "FEMesh.h"
#include "FEModelParam.h"

DEFINE_META_CLASS(FEBodyLoad, FEModelLoad, "");

//-----------------------------------------------------------------------------
FEBodyLoad::FEBodyLoad() : FEModelLoad()
{
}

//-----------------------------------------------------------------------------
FEBodyLoad::~FEBodyLoad()
{
}

//-----------------------------------------------------------------------------
//! initialization
bool FEBodyLoad::Init()
{
	// If the domain list is empty, add all the domains
	if (m_dom.Size() == 0)
	{
		FEMesh& mesh = GetMesh();
		for (int i=0; i<mesh.Domains(); ++i)
		{
			RgDomain* dom = &mesh.Domain(i);
			m_dom.Add(dom);
		}
	}
	return FEModelLoad::Init();
}

//-----------------------------------------------------------------------------
int FEBodyLoad::Domains() const
{
	return m_dom.Size();
}

//-----------------------------------------------------------------------------
RgDomain* FEBodyLoad::Domain(int i)
{
	return m_dom.Get(i);
}

//-----------------------------------------------------------------------------
void FEBodyLoad::SetDomainList(RgElementSet* elset)
{
	m_dom = *elset->GetDomainList();

	// add it to all the mapped parameters
	FEParameterList& PL = GetParameterList();
	FEParamIterator it = PL.first();
	for (int i = 0; i < PL.Parameters(); ++i, ++it)
	{
		FEParam& pi = *it;
		if (pi.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& param = pi.value<FEParamDouble>();
			param.SetItemList(elset);
		}
		else if (pi.type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& param = pi.value<FEParamVec3>();
			param.SetItemList(elset);
		}
		else if (pi.type() == FE_PARAM_MAT3D_MAPPED)
		{
			FEParamMat3d& param = pi.value<FEParamMat3d>();
			param.SetItemList(elset);
		}
	}
}

// get the domain list
RgDomainList& FEBodyLoad::GetDomainList()
{ 
	return m_dom; 
}

//! Serialization
void FEBodyLoad::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	m_dom.Serialize(ar);
}
