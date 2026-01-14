#include "RgTrussDomain.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "RgAssembler.h"

BEGIN_PARAM_DEFINE(RgTrussDomain, RgDomain)
	ADD_PARAMETER(m_a0, "cross_sectional_area");
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
//! Constructor
RgTrussDomain::RgTrussDomain(FEModel* pfem) : RgDomain(FE_DOMAIN_TRUSS, pfem), FEElasticDomain(pfem), m_dofU(pfem)
{
	m_a0 = 0.0;
	m_pMat = 0;
	m_assembler = nullptr;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		//m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
//! copy operator
RgTrussDomain& RgTrussDomain::operator = (RgTrussDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
bool RgTrussDomain::Create(int nsize, FE_Element_Spec espec)
{
	m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i) m_Elem[i].SetDomain(this);
	return true;
}

//-----------------------------------------------------------------------------
void RgTrussDomain::ForEachElement(std::function<void(RgElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}

//-----------------------------------------------------------------------------
void RgTrussDomain::ForEachTrussElement(std::function<void(FETrussElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(m_Elem[i]);
}

//-----------------------------------------------------------------------------
//! get the dof list
const FEDofList& RgTrussDomain::GetDOFList() const
{
	return m_dofU;
}

//-----------------------------------------------------------------------------
void RgTrussDomain::SetMaterial(FEMaterial* pmat)
{
	RgDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FETrussMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! initialize the domain
bool RgTrussDomain::Init()
{
	// base class initialization
	if (RgDomain::Init() == false) return false;

	// initialize truss data
	if (m_a0 != 0.0)
	{
		for (int i = 0; i < Elements(); ++i)
		{
			FETrussElement& el = Element(i);
			el.m_a0 = m_a0;
		}
	}

	for (int i = 0; i < (int)m_Elem.size(); ++i)
	{
		// unpack the element
		FETrussElement& el = m_Elem[i];

		// nodal coordinates
		Vector3d r0[2];
		for (int j = 0; j < 2; ++j) r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;

		// initial length
		el.m_L0 = (r0[1] - r0[0]).norm();
	}

	return true;
}

//-----------------------------------------------------------------------------
void RgTrussDomain::Reset()
{
	ForEachMaterialPoint([](RgMaterialPoint& mp) {
		mp.Init();
	});
}

//-----------------------------------------------------------------------------
void RgTrussDomain::UnpackLM(RgElement &el, vector<int>& lm)
{
	lm.resize(6);
	FENode& n1 = m_pMesh->Node(el.m_node[0]);
	FENode& n2 = m_pMesh->Node(el.m_node[1]);
	lm[0] = n1.m_dofs[m_dofU[0]];
	lm[1] = n1.m_dofs[m_dofU[1]];
	lm[2] = n1.m_dofs[m_dofU[2]];
	lm[3] = n2.m_dofs[m_dofU[0]];
	lm[4] = n2.m_dofs[m_dofU[1]];
	lm[5] = n2.m_dofs[m_dofU[2]];
}

//-----------------------------------------------------------------------------
void RgTrussDomain::Activate()
{
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void RgTrussDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.m_ra = m_pMesh->Node(el.m_node[0]).m_rt;
		el.m_rb = m_pMesh->Node(el.m_node[1]).m_rt;
	}
}

//-----------------------------------------------------------------------------
void RgTrussDomain::Update(const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
void RgTrussDomain::InternalForces(FEGlobalVector& R)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InternalForces(R);
		return;
	}
	
	// fall back to direct implementation
	vector<int> lm;

	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];

		// get the element force vector
		vector<double> fe;
		ElementInternalForces(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void RgTrussDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->StiffnessMatrix(LS);
		return;
	}
	
	// fall back to direct implementation
	vector<int> lm;
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];

		// element stiffness matrix
		matrix ke;
		ElementStiffness(i, ke);

		// get the element's LM vector
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void RgTrussDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->MassMatrix(LS, scale);
		return;
	}
	
	// fall back to direct implementation
	vector<int> lm;
	for (int i = 0; i < (int)m_Elem.size(); ++i)
	{
		FETrussElement& el = m_Elem[i];

		// element stiffness matrix
		matrix ke;
		ElementMassMatrix(el, ke);

		// get the element's LM vector
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates the truss element stiffness matrix
void RgTrussDomain::ElementStiffness(int iel, matrix& ke)
{
	FETrussElement& el = m_Elem[iel];

	// get the material
	FETrussMaterial* pmat = dynamic_cast<FETrussMaterial*>(GetMaterial());
	double E = pmat->m_E;

	// get the nodes
	vec3d& r1 = el.m_ra;
	vec3d& r2 = el.m_rb;
	double L = (r2 - r1).norm();

	// get the cross-sectional area
	double A = el.m_a0;

	// calculate the stiffness component
	double k = E*A/L;

	ke.resize(6, 6);
	ke.zero();

	// --- S O L I D ---
	double c12 = k / (L*L);
	double cx = (r2.x - r1.x);
	double cy = (r2.y - r1.y);
	double cz = (r2.z - r1.z);

	ke[0][0] = c12*cx*cx; ke[0][1] = c12*cx*cy; ke[0][2] = c12*cx*cz;
	ke[1][0] = c12*cy*cx; ke[1][1] = c12*cy*cy; ke[1][2] = c12*cy*cz;
	ke[2][0] = c12*cz*cx; ke[2][1] = c12*cz*cy; ke[2][2] = c12*cz*cz;

	ke[3][3] = c12*cx*cx; ke[3][4] = c12*cx*cy; ke[3][5] = c12*cx*cz;
	ke[4][3] = c12*cy*cx; ke[4][4] = c12*cy*cy; ke[4][5] = c12*cy*cz;
	ke[5][3] = c12*cz*cx; ke[5][4] = c12*cz*cy; ke[5][5] = c12*cz*cz;

	ke[0][3] = -ke[0][0]; ke[0][4] = -ke[0][1]; ke[0][5] = -ke[0][2];
	ke[1][3] = -ke[1][0]; ke[1][4] = -ke[1][1]; ke[1][5] = -ke[1][2];
	ke[2][3] = -ke[2][0]; ke[2][4] = -ke[2][1]; ke[2][5] = -ke[2][2];

	ke[3][0] = -ke[0][0]; ke[3][1] = -ke[0][1]; ke[3][2] = -ke[0][2];
	ke[4][0] = -ke[1][0]; ke[4][1] = -ke[1][1]; ke[4][2] = -ke[1][2];
	ke[5][0] = -ke[2][0]; ke[5][1] = -ke[2][1]; ke[5][2] = -ke[2][2];
}

//-----------------------------------------------------------------------------
//! Calculates the internal stress vector for solid elements
void RgTrussDomain::ElementInternalForces(FETrussElement& el, vector<double>& fe)
{
	// get the material
	FETrussMaterial* pmat = dynamic_cast<FETrussMaterial*>(GetMaterial());

	// get the nodes
	vec3d& r1 = el.m_ra;
	vec3d& r2 = el.m_rb;
	vec3d e12 = r2 - r1;
	double L = e12.norm();

	// get the cross-sectional area
	double A = el.m_a0;

	// calculate the force component
	double k = pmat->m_E*A/L;
	double f = k*(L - el.m_L0);

	// get the direction
	vec3d q12 = e12/L;

	// calculate internal force vector
	fe.resize(6);

	fe[0] =  f*q12.x;
	fe[1] =  f*q12.y;
	fe[2] =  f*q12.z;

	fe[3] = -f*q12.x;
	fe[4] = -f*q12.y;
	fe[5] = -f*q12.z;
}

//-----------------------------------------------------------------------------
void RgTrussDomain::ElementMassMatrix(FETrussElement& el, Matrix& ke)
{
	// get the material
	FETrussMaterial* pmat = dynamic_cast<FETrussMaterial*>(GetMaterial());

	// get the nodes
	vec3d& r1 = el.m_ra;
	vec3d& r2 = el.m_rb;
	double L = (r2 - r1).norm();

	// get the cross-sectional area
	double A = el.m_a0;

	// calculate the stiffness component
	double k = pmat->Density()*A*L/6.0;

	ke.resize(6, 6);
	ke.zero();

	ke[0][0] = 2*k; ke[0][1] = 0; ke[0][2] = 0;
	ke[1][0] = 0; ke[1][1] = 2*k; ke[1][2] = 0;
	ke[2][0] = 0; ke[2][1] = 0; ke[2][2] = 2*k;

	ke[3][3] = 2*k; ke[3][4] = 0; ke[3][5] = 0;
	ke[4][3] = 0; ke[4][4] = 2*k; ke[4][5] = 0;
	ke[5][3] = 0; ke[5][4] = 0; ke[5][5] = 2*k;

	ke[0][3] = k; ke[0][4] = 0; ke[0][5] = 0;
	ke[1][3] = 0; ke[1][4] = k; ke[1][5] = 0;
	ke[2][3] = 0; ke[2][4] = 0; ke[2][5] = k;

	ke[3][0] = k; ke[3][1] = 0; ke[3][2] = 0;
	ke[4][0] = 0; ke[4][1] = k; ke[4][2] = 0;
	ke[5][0] = 0; ke[5][1] = 0; ke[5][2] = k;
}

//-----------------------------------------------------------------------------
Vector3d RgTrussDomain::TrussNormal(FETrussElement& el)
{
	vec3d r1 = m_pMesh->Node(el.m_node[0]).m_rt;
	vec3d r2 = m_pMesh->Node(el.m_node[1]).m_rt;

	vec3d a = r2 - r1;
	double L = a.unit();

	// get the cross-sectional area
	double A = el.m_a0;

	// normal force
	FETrussMaterial* pmat = dynamic_cast<FETrussMaterial*>(GetMaterial());
	double F = pmat->m_E*A*(L - el.m_L0) / L;

	return a*(F / L);
}