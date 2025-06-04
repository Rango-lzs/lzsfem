#include "FEPrescribedDOF.h"
#include "FENodeSet.h"
#include "basicio/DumpStream.h"
#include "FEMesh.h"
#include "logger/log.h"
#include "materials/FEMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEPrescribedDOF, FEPrescribedNodeSet)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_dof  , "dof", 0, "$(dof_list)");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF() : FEPrescribedNodeSet()
{
	m_scale = 0.0;
	m_dof = -1;
}

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF(FEModel* pfem, int dof, FENodeSet* nset) : FEPrescribedNodeSet()
{
	m_scale = 0.0;
	SetNodeSet(nset);
	SetDOF(dof);
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::SetDOF(int ndof)
{
	m_dof = ndof;
}

//-----------------------------------------------------------------------------
bool FEPrescribedDOF::SetDOF(const char* szdof)
{
	int ndof = GetDOFIndex(szdof);
	assert(ndof >= 0);
	if (ndof < 0) return false;
	SetDOF(ndof);
	return true;
}

//-----------------------------------------------------------------------------
// Sets the displacement scale factor. An optional load curve index can be given
// of the load curve that will control the scale factor.
FEPrescribedDOF& FEPrescribedDOF::SetScale(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		AttachLoadController(&m_scale, lc);
	}
	return *this;
}

//-----------------------------------------------------------------------------
bool FEPrescribedDOF::Init()
{
	// set the dof first before calling base class
	if (m_dof < 0) return false;
	SetDOFList(m_dof);

	// don't forget to call the base class
	if (FEPrescribedNodeSet::Init() == false) return false;

	// make sure this is not a rigid node
	FEMesh& mesh = GetMesh();
	int NN = mesh.Nodes();
	const FENodeSet& nset = *GetNodeSet();
	for (size_t i = 0; i<nset.Size(); ++i)
	{
		int nid = nset[i];
		if ((nid < 0) || (nid >= NN)) return false;
		if (mesh.Node(nid).m_rid != -1)
		{
			feLogError("Rigid nodes cannot be prescribed.");
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::GetNodalValues(int n, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[n];
	const FENode& node = *nset.Node(n);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = n;

	val[0] = m_scale(mp);
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::CopyFrom(FEBoundaryCondition* pbc)
{
	FEPrescribedDOF* ps = dynamic_cast<FEPrescribedDOF*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	//CopyParameterListState(ps->GetParameterList());
}
