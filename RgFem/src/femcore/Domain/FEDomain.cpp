#include "femcore/FEDomain.h"
#include "materials/FEMaterial.h"
#include "basicio/DumpStream.h"
#include "femcore/FEMesh.h"
#include "femcore/Matrix/FEGlobalMatrix.h"

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int nclass, FEModel* fem) : FEMeshPartition(nclass, fem)
{

}

//-----------------------------------------------------------------------------
void FEDomain::SetMaterial(FEMaterial* pm)
{
	assert(pm);
	if (pm) pm->AddDomain(this);
}

//-----------------------------------------------------------------------------
void FEDomain::SetMatID(int mid)
{
	ForEachElement([=](FEElement& el) { el.SetMatID(mid); });
}

//-----------------------------------------------------------------------------
// This routine allocates the material point data for the element's integration points.
// Currently, this has to be called after the elements have been assigned a type (since this
// determines how many integration points an element gets). 
void FEDomain::CreateMaterialPointData()
{
	FEMaterial* pmat = GetMaterial();
	FEMesh* mesh = GetMesh();
	if (pmat) ForEachElement([=](FEElement& el) {

		Vector3d r[FEElement::MAX_NODES];
		int ne = el.Nodes();
		for (int i = 0; i < ne; ++i) r[i] = mesh->Node(el.m_node[i]).m_r0;

		for (int k = 0; k < el.GaussPoints(); ++k)
		{
			FEMaterialPoint* mp = new FEMaterialPoint(pmat->CreateMaterialPointData());
			mp->m_r0 = el.Evaluate(r, k);
			mp->m_index = k;
			el.SetMaterialPointData(mp, k);
		}
	});
}

//-----------------------------------------------------------------------------
// serialization
void FEDomain::Serialize(DumpStream& ar)
{
	FEMeshPartition::Serialize(ar);

	if (ar.IsShallow())
	{
		int NEL = Elements();
		for (int i = 0; i < NEL; ++i)
		{
			FEElement& el = ElementRef(i);
			el.Serialize(ar);
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			FEMaterial* mat = GetMaterial();
			ar << mat;

			int NEL = Elements();
			ar << NEL;
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			FEMaterial* pmat = 0;
			ar >> pmat;
			SetMaterial(pmat);

			FE_Element_Spec espec; // invalid element spec!

			int NEL = 0;
			ar >> NEL;
			Create(NEL, espec);
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j)
				{
					FEMaterialPoint* mp = new FEMaterialPoint(pmat->CreateMaterialPointData());
					el.SetMaterialPointData(mp, j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Unpack the LM data for an element of this domain
void FEDomain::UnpackLM(FEElement& el, std::vector<int>& lm)
{
	UnpackLM(el, GetDOFList(), lm);
}

//-----------------------------------------------------------------------------
//! Activate the domain
void FEDomain::Activate()
{
	Activate(GetDOFList());
}

//-----------------------------------------------------------------------------
// This is the default packing method. 
// It stores all the degrees of freedom for the first node in the order defined
// by the DOF array, then for the second node, and so on. 
void FEDomain::UnpackLM(FEElement& ele, const FEDofList& dof, std::vector<int>& lm)
{
	FEMesh* mesh = GetMesh();
	int N = ele.Nodes();
	int ndofs = dof.Size();
	lm.resize(N*ndofs);
	for (int i = 0; i<N; ++i)
	{
		FENode& node = ele.giveNode(i);
		std::vector<int>& id = node.getDofs();
		for (int j = 0; j<ndofs; ++j) lm[i*ndofs + j] = id[dof[j]];
	}
}

//-----------------------------------------------------------------------------
void FEDomain::BuildMatrixProfile(FEGlobalMatrix& M)
{
	std::vector<int> elm;
	const int NE = Elements();
	for (int j = 0; j<NE; ++j)
	{
		FEElement& el = ElementRef(j);
		UnpackLM(el, elm);
		M.build_add(elm);
	}
}

//-----------------------------------------------------------------------------
void FEDomain::Activate(const FEDofList& dof)
{
	// get the number of degrees of freedom for this domain.
	const int ndofs = dof.Size();

	// activate all the degrees of freedom of this domain
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			for (int j = 0; j<ndofs; ++j)
				if (dof[j] >= 0) node.set_active(dof[j]);
		}
	}
}
