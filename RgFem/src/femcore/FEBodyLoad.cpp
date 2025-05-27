
#include "FEBodyLoad.h"
#include "FEMesh.h"
#include "FEModelParam.h"

DEFINE_META_CLASS(FEBodyLoad, FEModelLoad, "");

//-----------------------------------------------------------------------------
FEBodyLoad::FEBodyLoad(FEModel* pfem) : FEModelLoad(pfem)
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
	if (m_dom.IsEmpty())
	{
		FEMesh& mesh = GetMesh();
		for (int i=0; i<mesh.Domains(); ++i)
		{
			FEDomain* dom = &mesh.Domain(i);
			m_dom.AddDomain(dom);
		}
	}
	return FEModelLoad::Init();
}

//-----------------------------------------------------------------------------
int FEBodyLoad::Domains() const
{
	return m_dom.Domains();
}

//-----------------------------------------------------------------------------
FEDomain* FEBodyLoad::Domain(int i)
{
	return m_dom.GetDomain(i);
}

//-----------------------------------------------------------------------------
void FEBodyLoad::SetDomainList(FEElementSet* elset)
{
	m_dom = elset->GetDomainList();

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
FEDomainList& FEBodyLoad::GetDomainList()
{ 
	return m_dom; 
}

//! Serialization
void FEBodyLoad::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	m_dom.Serialize(ar);
}
