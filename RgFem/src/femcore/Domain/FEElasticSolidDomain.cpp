#include "FEElasticSolidDomain.h"
#include "materials/FEElasticMaterial.h"
#include "femcore/FEBodyForce.h"
#include "logger/log.h"
#include <femcore/FEModel.h>
#include "femcore/FEAnalysis/FEAnalysis.h"
#include "femcore/sys.h"
#include "femcore/FELinearSystem.h"
#include "materials/FEElasticMaterialPoint.h"
#include "../FEException.h"
#include "../FEMeshPartition.h"
#include "femcore/FEMesh.h"
#include "../FESolidModule.h"

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEElasticSolidDomain::FEElasticSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dofSU(pfem), m_dofV(pfem), m_dofSV(pfem), m_dofSA(pfem), m_dof(pfem)
{
	m_pMat = 0;
    m_alphaf = m_beta = 1;
    m_alpham = 2;
	m_update_dynamic = true; // default for backward compatibility

	m_secant_stress = false;
	m_secant_tangent = false;

	// TODO: Can this be done in Init, since  there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(GetVariableName(DISPLACEMENT));
		m_dofR.AddVariable(GetVariableName(RIGID_ROTATION));
		m_dofSU.AddVariable(GetVariableName(SHELL_DISPLACEMENT));
		m_dofV.AddVariable(GetVariableName(VELOCTIY));
		m_dofSV.AddVariable(GetVariableName(SHELL_VELOCITY));
		m_dofSA.AddVariable(GetVariableName(SHELL_ACCELERATION));
	}
}

//-----------------------------------------------------------------------------
// get the total dof list
const FEDofList& FEElasticSolidDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEElasticSolidDomain& FEElasticSolidDomain::operator = (FEElasticSolidDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
//! Set flag for update for dynamic quantities
void FEElasticSolidDomain::SetDynamicUpdateFlag(bool b)
{
	m_update_dynamic = b;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEElasticSolidDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	if (pmat)
	{
		m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
		assert(m_pMat);
	}
	else m_pMat = 0;
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
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
//! serialization
void FEElasticSolidDomain::Serialize(DumpStream& ar)
{
	//erialize the base class, which instantiates the elements
	FESolidDomain::Serialize(ar);
	if (ar.IsShallow()) return;

	// serialize class variables
	ar & m_alphaf;
	ar & m_alpham;
	ar & m_beta;
	ar & m_update_dynamic;
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEElasticSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    m_alphaf = timeInfo.alphaf;
    m_alpham = timeInfo.alpham;
    m_beta = timeInfo.beta;

#pragma omp parallel for
	for (int i=0; i<Elements(); ++i)
	{
		FESolidElement& el = *m_Elem[i];
		if (el.isActive())
		{
			int n = el.GaussPointSize();
			for (int j = 0; j < n; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				pt.m_Wp = pt.m_Wt;

				mp.Update(timeInfo);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::InternalForces(FEGlobalVector& R)
{
	int NE = Elements();
	#pragma omp parallel for shared (NE)
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = *m_Elem[i];

		if (1/*el.isActive()*/) {  //Element需要有isActive的标识吗？
			// element force std::vector
			std::vector<double> fe;
			std::vector<int> lm;

			// get the element force std::vector and initialize it to zero
			int ndof = 3 * el.NodeSize();
			fe.assign(ndof, 0);

			// calculate internal force std::vector
			ElementInternalForce(el, fe);

			// get the element's LM std::vector
			UnpackLM(el, lm);

			// assemble element 'fe'-std::vector into global R std::vector
			// Rango TODO:
			//R.Assemble(el.getNodeIds(), lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
//calculates the internal equivalent nodal forces for solid elements
//B*sigma
void FEElasticSolidDomain::ElementInternalForce(FESolidElement& el, std::vector<double>& fe)
{
	// jacobian Matrix, inverse jacobian Matrix and determinants
	double Ji[3][3];

	int nint = el.GaussPointSize();
	int neln = el.NodeSize();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		double detJt = (m_update_dynamic ? invjact(el, Ji, n, m_alphaf) : invjact(el, Ji, n));

		detJt *= gw[n];

		// get the stress std::vector for this integration point
        const Matrix3ds& s = pt.m_s;  //Cauchy stress

		const double* Gr = el.Gr(n);
		const double* Gs = el.Gs(n);
		const double* Gt = el.Gt(n);

		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual std::vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
	// define some parameters that will be passed to lambda
	FEBodyForce* bodyForce = &BF;

	// TODO: a remaining issue here is that dofU does not consider the shell displacement
	// dofs for interface nodes (see UnpackLM). Is that an issue?

	// evaluate the residual contribution
	LoadVector(R, m_dofU, [=](FEMaterialPoint& mp, int node_a, std::vector<double>& fa) {

		// evaluate density
		double density = m_pMat->Density(mp);

		// get the force
		Vector3d f = bodyForce->force(mp);

		// get element shape functions
		double* H = mp.m_shape;

		// get the initial Jacobian
		double J0 = mp.m_J0;

		// set integrand
		fa[0] = -H[node_a] * density* f.x * J0;
		fa[1] = -H[node_a] * density* f.y * J0;
		fa[2] = -H[node_a] * density* f.z * J0;
	});
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n
void FEElasticSolidDomain::ElementGeometricalStiffness(FESolidElement &el, Matrix &ke)
{
	// spatial derivatives of shape functions
	Vector3d G[FEElement::MAX_NODES];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate geometrical element stiffness Matrix
	int neln = el.NodeSize();
	int nint = el.GaussPointSize();
	for (int n = 0; n<nint; ++n)
	{
		// calculate shape function gradients and jacobian
		double w = ShapeGradient(el, n, G, m_alphaf)*gw[n]*m_alphaf;

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		Matrix3ds& s = pt.m_s;

		for (int i = 0; i<neln; ++i)
			for (int j = 0; j<neln; ++j)
			{
				double kab = (G[i]*(s * G[j]))*w;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element Matrix

void FEElasticSolidDomain::ElementMaterialStiffness(FESolidElement &el, Matrix &ke)
{
	// Get the current element's data
	const int nint = el.GaussPointSize();
	const int neln = el.NodeSize();

	// global derivatives of shape functions
	Vector3d G[FEElement::MAX_NODES];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' Matrix
	double D[6][6] = {0};	// The 'D' Matrix

	// The 'D*BL' Matrix
	double DBL[6][3];

	// jacobian
	double detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness Matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian and shape function gradients
        // alpha,HHT time integration parameters
		//G ： dH/dx
		detJt = ShapeGradient(el, n, G, m_alphaf)*gw[n]*m_alphaf;

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		// get the 'D' Matrix
//		tens4ds C = m_pMat->Tangent(mp);
        tens4dmm C = (m_secant_tangent ? m_pMat->SecantTangent(mp) : m_pMat->SolidTangent(mp));
		C.extract(D);

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		//按照3*3的分块矩阵进行计算
		for (int i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = G[i].x;
			Gyi = G[i].y;
			Gzi = G[i].z;

			for (int j=0, j3 = 0; j<neln; ++j, j3 += 3)
			{
				Gxj = G[j].x;
				Gyj = G[j].y;
				Gzj = G[j].z;

				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;

				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;

				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::StiffnessMatrix(FELinearSystem& ls)
{
	// repeat over all solid elements
	int NE = Elements();
	
	#pragma omp parallel for shared (NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& elem = *m_Elem[iel];

		if (elem.isActive()) {

			// get the element's LM std::vector
			std::vector<int> lm;
			UnpackLM(elem, lm);

			// element stiffness Matrix
			FEElementMatrix ke(elem, lm);

			// create the element's stiffness Matrix
			int ndof = 3 * elem.NodeSize();
			ke.resize(ndof, ndof);
			ke.zero();

			// calculate geometrical stiffness
			ElementGeometricalStiffness(elem, ke);

			// calculate material stiffness
			ElementMaterialStiffness(elem, ke);

/*			// assign symmetic parts
			// TODO: Can this be omitted by changing the Assemble routine so that it only
			// grabs elements from the upper diagonal Matrix?
			for (int i = 0; i < ndof; ++i)
				for (int j = i + 1; j < ndof; ++j)
					ke[j][i] = ke[i][j];
*/
			// assemble element Matrix in global stiffness Matrix
			ls.Assemble(ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// TODO: a remaining issue here is that dofU does not consider the shell displacement
	// dofs for interface nodes (see UnpackLM). Is that an issue?

	// evaluate body force stiffness
	LoadStiffness(LS, m_dofU, m_dofU, [=](FEMaterialPoint& mp, int node_a, int node_b, Matrix& Kab) {

		// density
		double density = m_pMat->Density(mp);

		// shape functions
		double* H = mp.m_shape;

		// Jacobian
		double J0 = mp.m_J0;

		// mass
		double kab = scale *density*H[node_a] * H[node_b] * J0;
		Kab.zero();
		Kab[0][0] = kab;
		Kab[1][1] = kab;
		Kab[2][2] = kab;
	});
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// define some parameters that will be passed to lambda
	FESolidMaterial* mat = m_pMat;
	FEBodyForce* bodyForce = &bf;

	// TODO: a remaining issue here is that dofU does not consider the shell displacement
	// dofs for interface nodes (see UnpackLM). Is that an issue?

	// evaluate body force stiffness
	LoadStiffness(LS, m_dofU, m_dofU, [=](FEMaterialPoint& mp, int node_a, int node_b, Matrix& Kab) {

		// loop over integration points
		double detJ = mp.m_J0 * m_alphaf;

		// density
		double dens_n = mat->Density(mp);
			
		// get the stiffness
		Matrix3d K = bodyForce->stiffness(mp);

		// shape functions
		double* H = mp.m_shape;

		// put it together
		Kab = K*(-H[node_a] * H[node_b] * dens_n*detJ);
	});
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness Matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal Matrix due to the symmetry of the element stiffness Matrix
//! The last section of this function fills the rest of the element stiffness Matrix.

void FEElasticSolidDomain::ElementStiffness(const FETimeInfo& tp, int iel, Matrix& ke)
{
	FESolidElement& el = Element(iel);

	// calculate material stiffness (i.e. constitutive component)
	ElementMaterialStiffness(el, ke);

	// calculate geometrical stiffness
	ElementGeometricalStiffness(el, ke);

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal Matrix?
	int ndof = 3*el.NodeSize();
	int i, j;
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Update(const FETimeInfo& tp)
{
	bool berr = false;
	int NE = Elements();
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			FESolidElement& el = Element(i);
			if (el.isActive())
			{
				UpdateElementStress(i, tp);
			}
		}
		catch (NegativeJacobian e)
		{
			#pragma omp critical
			{
				// reset the logfile mode
				berr = true;
				if (e.DoOutput()) feLogError(e.what());
			}
		}
	}

	if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FEElasticSolidDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double dt =tp.timeIncrement;
    
	// get the solid element
	FESolidElement& el = *m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPointSize();

	// number of nodes
	int neln = el.NodeSize();

	// nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    Vector3d r[NELN], v[NELN], a[NELN];
	GetCurrentNodalCoordinates(el, r, m_alphaf);

	// update dynamic quantities
	if (m_update_dynamic)
	{
		for (int j = 0; j<neln; ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			v[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*m_alphaf + node.m_vp*(1 - m_alphaf);
			a[j] = node.m_at*m_alpham + node.m_ap*(1 - m_alpham);
		}
	}

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// material point coordinates
		mp.m_rt = el.Evaluate(r, n);

		// get the deformation gradient and determinant at intermediate time
        double Jt;
        Matrix3d Ft, Fp;
        Jt = defgrad(el, Ft, n);
        defgradp(el, Fp, n);

		if (m_alphaf == 1.0)
		{
			pt.m_F = Ft;
            pt.m_J = Jt;
		}
		else
		{
			pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
            pt.m_J = pt.m_F.det();
		}

        Matrix3d Fi = pt.m_F.inverse();
        pt.m_L = (Ft - Fp)*Fi / dt;
		if (m_update_dynamic)
		{
			pt.m_v = el.Evaluate(v, n);
			pt.m_a = el.Evaluate(a, n);
		}

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, tp);
        
		// calculate the stress at this material point
//		pt.m_s = m_pMat->Stress(mp);
		pt.m_s = (m_secant_stress ? m_pMat->SecantStress(mp) : m_pMat->Stress(mp));
        
        // adjust stress for strain energy conservation
        if (m_alphaf == 0.5) 
		{
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(m_pMat);

			// evaluate strain energy at current time
			Matrix3d Ftmp = pt.m_F;
			double Jtmp = pt.m_J;
			pt.m_F = Ft;
			pt.m_J = Jt;
			pt.m_Wt = pme->StrainEnergyDensity(mp);
			pt.m_F = Ftmp;
			pt.m_J = Jtmp;

            Matrix3ds D = pt.RateOfDeformation();
            double D2 = D.dotdot(D);
            if (D2 > 0)
                pt.m_s += D*(((pt.m_Wt-pt.m_Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEElasticSolidDomain::UnpackLM(FEElement& el, std::vector<int>& lm)
{
	int N = el.NodeSize();
	lm.resize(N*6);
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		std::vector<int>& id = node.m_dofs;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofU[0]];
		lm[3*i+1] = id[m_dofU[1]];
		lm[3*i+2] = id[m_dofU[2]];

		// rigid rotational dofs
		lm[3*N + 3*i  ] = id[m_dofR[0]];
		lm[3*N + 3*i+1] = id[m_dofR[1]];
		lm[3*N + 3*i+2] = id[m_dofR[2]];
	}
    
    // substitute interface dofs for solid-shell interfaces
	FESolidElement& sel = static_cast<FESolidElement&>(el);
	for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            std::vector<int>& id = node.m_dofs;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofSU[0]];
            lm[3*i+1] = id[m_dofSU[1]];
            lm[3*i+2] = id[m_dofSU[2]];
        }
    }
}

//-----------------------------------------------------------------------------
// Calculate inertial forces \todo Why is F no longer needed?
void FEElasticSolidDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
    int NE = Elements();
#pragma omp parallel for shared(R, F)
	for (int i=0; i<NE; ++i)
    {
		// get the element
		FESolidElement& el = *m_Elem[i];

		if (el.isActive()) {
			// element force std::vector
			std::vector<double> fe;
			std::vector<int> lm;

			// get the element force std::vector and initialize it to zero
			int ndof = 3 * el.NodeSize();
			fe.assign(ndof, 0);

			// calculate internal force std::vector
			ElementInertialForce(el, fe);

			// get the element's LM std::vector
			UnpackLM(el, lm);

			// assemble element 'fe'-std::vector into global R std::vector
			R.Assemble(el.m_node, lm, fe);
		}
    }
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::ElementInertialForce(FESolidElement& el, std::vector<double>& fe)
{
    int nint = el.GaussPointSize();
    int neln = el.NodeSize();
    
    double*    gw = el.GaussWeights();

    // repeat for all integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        double dens = m_pMat->Density(mp);
        double J0 = detJ0(el, n)*gw[n];
        
        double* H = el.H(n);
        for (int i=0; i<neln; ++i)
        {
            double tmp = H[i]*J0*dens;
            fe[3*i  ] -= tmp*pt.m_a.x;
            fe[3*i+1] -= tmp*pt.m_a.y;
            fe[3*i+2] -= tmp*pt.m_a.z;
        }
    }
}


//=================================================================================================

BEGIN_PARAM_DEFINE(FEStandardElasticSolidDomain, FEElasticSolidDomain)
	ADD_PARAMETER(m_elemType, "elem_type", FE_PARAM_ATTRIBUTE, "$(solid_element)\0");
END_PARAM_DEFINE();

FEStandardElasticSolidDomain::FEStandardElasticSolidDomain(FEModel* fem) : FEElasticSolidDomain(fem)
{

}
