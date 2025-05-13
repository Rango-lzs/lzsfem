#include "femcore/Domain/FEShellDomain.h"
#include "femcore/FEMesh.h"
#include "materials/FEMaterial.h"

//-----------------------------------------------------------------------------
//! constructor
FEShellDomain::FEShellDomain(FEModel* fem) : FEDomain(FE_DOMAIN_SHELL, fem)
{
}

//-----------------------------------------------------------------------------
void FEShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	ForEachShellElement([](FEShellElement& el) {
		int ni = el.GaussPointSize();
		for (int j = 0; j<ni; ++j) el.GetMaterialPoint(j)->Init();

		int ne = el.NodeSize();
		for (int j = 0; j<ne; ++j) el.m_ht[j] = el.m_h0[j];
	});
}

//-----------------------------------------------------------------------------
void FEShellDomain::InitShells()
{
	ForEachShellElement([](FEShellElement& el) {
		int n = el.NodeSize();
		for (int j = 0; j<n; ++j) el.m_ht[j] = el.m_h0[j];
	});
}

//-----------------------------------------------------------------------------
//! get the current nodal coordinates
void FEShellDomain::GetCurrentNodalCoordinates(const FEShellElement& el, Vector3d* rt, const bool back)
{
    int neln = el.NodeSize();
    if (!back)
        for (int i = 0; i<neln; ++i) rt[i] = m_pMesh->Node(el.getNodeId(i)).m_rt;
    else
        for (int i = 0; i<neln; ++i) rt[i] = m_pMesh->Node(el.getNodeId(i)).st();
}

//-----------------------------------------------------------------------------
//! get the current nodal coordinates
void FEShellDomain::GetCurrentNodalCoordinates(const FEShellElement& el, Vector3d* rt, double alpha, const bool back)
{
    int neln = el.NodeSize();
    if (!back) {
        for (int i = 0; i<neln; ++i) {
            FENode& nd = m_pMesh->Node(el.getNodeId(i));
            rt[i] = nd.m_rt*alpha + nd.m_rp*(1 - alpha);
        }
    }
    else {
        for (int i = 0; i<neln; ++i) {
            FENode& nd = m_pMesh->Node(el.getNodeId(i));
            rt[i] = nd.st()*alpha + nd.sp()*(1 - alpha);
        }
    }
}

//-----------------------------------------------------------------------------
//! get the reference nodal coordinates
void FEShellDomain::GetReferenceNodalCoordinates(const FEShellElement& el, Vector3d* r0, const bool back)
{
    int neln = el.NodeSize();
    if (!back)
        for (int i = 0; i<neln; ++i) r0[i] = m_pMesh->Node(el.getNodeId(i)).m_r0;
    else
        for (int i = 0; i<neln; ++i) r0[i] = m_pMesh->Node(el.getNodeId(i)).s0();
}

//-----------------------------------------------------------------------------
//! get the previous nodal coordinates
void FEShellDomain::GetPreviousNodalCoordinates(const FEShellElement& el, Vector3d* rp, const bool back)
{
    int neln = el.NodeSize();
    if (!back)
        for (int i = 0; i<neln; ++i) rp[i] = m_pMesh->Node(el.getNodeId(i)).m_rp;
    else
        for (int i = 0; i<neln; ++i) rp[i] = m_pMesh->Node(el.getNodeId(i)).sp();
}

//-----------------------------------------------------------------------------
void FEShellDomain::ForEachShellElement(std::function<void(FEShellElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(Element(i));
}

//=================================================================================================

FEShellDomainOld::FEShellDomainOld(FEModel* fem) : FEShellDomain(fem)
{
}

//-----------------------------------------------------------------------------
bool FEShellDomainOld::Create(int nelems, FE_Element_Spec espec)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementOld& el = m_Elem[i];
		el.setLocalId(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i=0; i<nelems; ++i) m_Elem[i].setType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
double FEShellDomainOld::Volume(FEShellElement& se)
{
	FEShellElementOld& el = static_cast<FEShellElementOld&>(se);

	int neln = el.NodeSize();

	// initial nodal coordinates and directors
	Vector3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (int i = 0; i<neln; ++i)
	{
		r0[i] = Node(el.getLocNodeId(i)).m_r0;
		D0[i] = el.m_D0[i];
	}

	int nint = el.GaussPointSize();
	double *w = el.GaussWeights();
	double V = 0;
	Vector3d g[3];
	for (int n = 0; n<nint; ++n)
	{
		// jacobian matrix
		double eta = el.gt(n);

		double* Mr = el.Hr(n);
		double* Ms = el.Hs(n);
		double* M = el.H(n);

		// evaluate covariant basis vectors
		g[0] = g[1] = g[2] = Vector3d(0, 0, 0);
		for (int i = 0; i<neln; ++i)
		{
			g[0] += (r0[i] + D0[i] * eta / 2)*Mr[i];
			g[1] += (r0[i] + D0[i] * eta / 2)*Ms[i];
			g[2] += D0[i] * (M[i] / 2);
		}

		Matrix3d J = Matrix3d(g[0].x, g[1].x, g[2].x,
			g[0].y, g[1].y, g[2].y,
			g[0].z, g[1].z, g[2].z);

		// calculate the determinant
		double detJ0 = J.det();

		V += detJ0*w[n];
	}

	return V;
}

//-----------------------------------------------------------------------------
//! Calculate all shell normals (i.e. the shell directors).
//! And find shell nodes
void FEShellDomainOld::InitShells()
{
	FEShellDomain::InitShells();

	FEMesh& mesh = *GetMesh();
	for (int i = 0; i<Elements(); ++i)
	{
		FEShellElementOld& el = ShellElement(i);
		int ne = el.NodeSize();
		for (int j = 0; j<ne; ++j)
		{
			Vector3d d0 = mesh.Node(el.getNodeId(j)).m_d0;
			d0.unit();
			el.m_D0[j] = d0 * el.m_h0[j];
		}
	}
}

//=================================================================================================

BEGIN_PARAM_DEFINE(FEShellDomainNew, FEShellDomain)
	ADD_PARAMETER(m_h0, "shell_thickness");
END_PARAM_DEFINE();

FEShellDomainNew::FEShellDomainNew(FEModel* fem) : FEShellDomain(fem)
{
	m_h0 = 0.0;
}

//-----------------------------------------------------------------------------
bool FEShellDomainNew::Create(int nelems, FE_Element_Spec espec)
{
	FEShellElementNew* p = new FEShellElementNew();
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementNew& el = m_Elem[i];
		el.setLocalId(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i = 0; i<nelems; ++i) m_Elem[i].setType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::AssignDefaultShellThickness()
{
	double h0 = DefaultShellThickness();
	if (h0 <= 0.0) return;

	for (int j = 0; j < Elements(); ++j)
	{
		FEShellElement& el = Element(j);
		int ne = el.NodeSize();
		for (int n = 0; n < ne; ++n) el.m_ht[n] = el.m_h0[n] = h0;
	}
}

//-----------------------------------------------------------------------------
double FEShellDomainNew::Volume(FEShellElement& se)
{
	FEShellElementNew& el = static_cast<FEShellElementNew&>(se);

	int neln = el.NodeSize();

	// initial nodal coordinates and directors
	Vector3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (int i = 0; i<neln; ++i)
	{
		r0[i] = Node(el.getLocNodeId(i)).m_r0;
		D0[i] = Node(el.getLocNodeId(i)).m_d0;
	}

	int nint = el.GaussPointSize();
	double *w = el.GaussWeights();
	double V = 0;
	Vector3d g[3];
	for (int n = 0; n<nint; ++n)
	{
		// jacobian matrix
		double eta = el.gt(n);

		double* Mr = el.Hr(n);
		double* Ms = el.Hs(n);
		double* M = el.H(n);

		// evaluate covariant basis vectors
		g[0] = g[1] = g[2] = Vector3d(0, 0, 0);
		for (int i = 0; i<neln; ++i)
		{
			g[0] += (r0[i] + D0[i] * eta / 2)*Mr[i];
			g[1] += (r0[i] + D0[i] * eta / 2)*Ms[i];
			g[2] += D0[i] * (M[i] / 2);
		}

		Matrix3d J = Matrix3d(g[0].x, g[1].x, g[2].x,
			g[0].y, g[1].y, g[2].y,
			g[0].z, g[1].z, g[2].z);

		// calculate the determinant
		double detJ0 = J.det();

		V += detJ0*w[n];
	}

	return V;
}


