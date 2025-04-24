#include "FETrussDomain.h"
#include "femcore/FEMesh.h"
#include "elements/FETrussElement.h"

//-----------------------------------------------------------------------------
FETrussDomain::FETrussDomain(FEModel* fem) : FEBeamDomain(fem)
{
}

//-----------------------------------------------------------------------------
bool FETrussDomain::Create(int nsize, FE_Element_Spec espec)
{
	m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.setLocalID(i);
		el.setMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i=0; i<nsize; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
Vector3d FETrussDomain::TrussNormal(FETrussElement& el)
{
	Vector3d rt[2];
	rt[0] = m_pMesh->Node(el.m_node[0]).m_rt;
	rt[1] = m_pMesh->Node(el.m_node[1]).m_rt;

	Vector3d a = rt[0];
	Vector3d b = rt[1];
	Vector3d n = b - a;
	n.unit();
	return n;
}

//-----------------------------------------------------------------------------
void FETrussDomain::ForEachTrussElement(std::function<void(FETrussElement& el)> f)
{
	int N = Elements();
	for (int i = 0; i < N; ++i)
	{
		FETrussElement& el = m_Elem[i];
		f(el);
	}
}
