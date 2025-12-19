#include "RgElasticSolidDomain.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "FEElasticSolidAssembler.h"

BEGIN_PARAM_DEFINE(RgElasticSolidDomain, RgSolidDomain)
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
//! Constructor
RgElasticSolidDomain::RgElasticSolidDomain(FEModel* pfem) : RgSolidDomain(pfem), RgElasticDomain(pfem), m_dofX(pfem)
{
	m_pMat = 0;
	m_assembler = 0;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		//m_dofX.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
//! copy operator
RgElasticSolidDomain& RgElasticSolidDomain::operator = (RgElasticSolidDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
bool RgElasticSolidDomain::Create(int nsize, FE_Element_Spec espec)
{
	m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i) m_Elem[i].SetDomain(this);
	return true;
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::ForEachElement(std::function<void(FEElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::ForEachSolidElement(std::function<void(FESolidElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(m_Elem[i]);
}

//-----------------------------------------------------------------------------
//! get the dof list
const FEDofList& RgElasticSolidDomain::GetDOFList() const
{
	return m_dofX;
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::SetMaterial(FEMaterial* pmat)
{
	RgSolidDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! initialize the domain
bool RgElasticSolidDomain::Init()
{
	// base class initialization
	if (RgSolidDomain::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::Reset()
{
	ForEachMaterialPoint([](FEMaterialPoint& mp) {
		mp.Init();
	});
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		lm[3*i  ] = node.m_dofs[m_dofX[0]];
		lm[3*i+1] = node.m_dofs[m_dofX[1]];
		lm[3*i+2] = node.m_dofs[m_dofX[2]];
	}
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofX[0]);
				node.set_active(m_dofX[1]);
				node.set_active(m_dofX[2]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::Update(const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::InternalForces(FEGlobalVector& R)
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
		FESolidElement& el = m_Elem[i];

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
void RgElasticSolidDomain::StiffnessMatrix(FELinearSystem& LS)
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
		FESolidElement& el = m_Elem[i];

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
void RgElasticSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
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
		FESolidElement& el = m_Elem[i];

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
//! calculates the solid element stiffness matrix
void RgElasticSolidDomain::ElementStiffness(int iel, matrix& ke)
{
	FESolidElement& el = m_Elem[iel];

	// Get the current element's data
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// get the material
	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial());
	assert(pme);

	// element stiffness matrix
	ke.resize(3*neln, 3*neln);
	ke.zero();

	// repeat for all integration points
	for (int i=0; i<nint; ++i)
	{
		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(i);
		double detJ = mp.determinant_of_jacobian;
		double w = mp.integration_weight;
		double Wp = detJ*w;

		// get the stiffness matrix
		tens4ds C = pme->Tangent(mp);

		// get the B-matrix
		Matrix3ds B;
		el.BlkStrainGradient(i, B);

		// multiply C:B
		Matrix3ds CB = C*B;

		// calculate element stiffness matrix
		for (int j=0; j<neln; ++j)
			for (int k=0; k<neln; ++k)
			{
				Matrix3d Kjk = CB*B(k,i);
				ke.add(3*j, 3*k, Kjk, Wp);
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the internal stress vector for solid elements
void RgElasticSolidDomain::ElementInternalForces(FESolidElement& el, vector<double>& fe)
{
	// get the material
	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial());
	assert(pme);

	// number of nodes
	int neln = el.Nodes();

	// number of integration points
	int nint = el.GaussPoints();

	// resize the force vector
	fe.resize(3*neln);

	// repeat for all integration points
	for (int i=0; i<nint; ++i)
	{
		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(i);
		double detJ = mp.determinant_of_jacobian;
		double w = mp.integration_weight;
		double Wp = detJ*w;

		// get the stress vector
		mat3ds s = pme->Stress(mp);

		// get the B-matrix
		Matrix3ds Bi;
		el.BlkStrainGradient(i, Bi);

		// calculate the internal force
		// fi = bi^T * si
		for (int j=0; j<neln; ++j)
		{
			vec3d fJ = Bi(j,i)*s;
			fe[3*j  ] -= fJ.x*Wp;
			fe[3*j+1] -= fJ.y*Wp;
			fe[3*j+2] -= fJ.z*Wp;
		}
	}
}

//-----------------------------------------------------------------------------
void RgElasticSolidDomain::ElementMassMatrix(FESolidElement& el, Matrix& ke)
{
	// Get the current element's data
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// element stiffness matrix
	ke.resize(3*neln, 3*neln);
	ke.zero();

	// repeat for all integration points
	for (int i=0; i<nint; ++i)
	{
		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(i);
		double detJ = mp.determinant_of_jacobian;
		double w = mp.integration_weight;
		double Wp = detJ*w;

		// get the density
		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial());
		double dens = pme->Density();

		// shape functions
		double* N = el.Gr(i);

		// calculate element stiffness matrix
		for (int j=0; j<neln; ++j)
			for (int k=0; k<neln; ++k)
			{
				double kab = N[j]*N[k]*dens*Wp;
				ke[3*j  ][3*k  ] += kab;
				ke[3*j+1][3*k+1] += kab;
				ke[3*j+2][3*k+2] += kab;
			}
	}
}