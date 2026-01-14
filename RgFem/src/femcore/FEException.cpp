#include "FEException.h"
#include "femcore/Solver/FESolver.h" // for FENodalDofInfo
#include <stdarg.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEException::FEException(const char* msg, int level)
{
	if (msg) m_what = msg;
	m_level = level;
}

FEException::~FEException()
{

}

const char* FEException::what()
{
	return m_what.c_str();
}

void FEException::what(const char* msg, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char sztxt[1024] = { 0 };
	va_start(args, msg);
	vsprintf(sztxt, msg, args);
	va_end(args);

	m_what = sztxt;
}

//-----------------------------------------------------------------------------
int FEException::level() const
{
	return m_level;
}

//-----------------------------------------------------------------------------
//! \todo implement error handling
ZeroDiagonal::ZeroDiagonal(int node, int ndof)
{
	what("Zero diagonal detected. Aborting run.");
}

/*
ZeroDiagonal::ZeroDiagonal(vector<int>& l, FEM& fem)
{
	// let's find what dof this equation belonges to
	FEMesh& m = fem.GetMesh();

	int nz = (int) l.size();
	int i, j, id;
	for (i=0; i<fem.m_nrb; ++i)
	{
		FERigidBody& rb = fem.m_RB[i];
		for (j=0; j<6; ++j)
		{
			id = rb.m_LM[j];
			if (id < -1) id = -id-2;

			for (int k=0; k<nz; ++k)
			{
				int n = l[k];
				if (id == n)
				{
					felog.printf("Zero diagonal on row %d.\nThis dof belongs to rigid body %d (dof %d)\n", n, i+1, j+1);
				}
			}
		}
	}

	vector<EQUATION> EQT(fem.m_neq);
	for (i=0; i<fem.m_neq; ++i)
	{
		EQUATION& q = EQT[i];
		q.node = -1;
		q.dof = -1;
	}

	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)
		{
			id = node.m_ID[j];
			if (id < -1) id = -id-2;
			if (id != -1)
			{
				assert(id < fem.m_neq);
				EQT[id].node = i;
				EQT[id].dof = j;
			}
		}
	}

	const int NMAX = 128;
	int N = (nz < NMAX? nz : NMAX), n;
	for (i=0; i<N; ++i)
	{
		n = l[i];
		EQUATION& q = EQT[n];
		felog.printf("Zero diagonal on row %d. (node %d, dof %d)\n", n+1, q.node+1, q.dof+1);
	}
	if (nz > NMAX) felog.printf("(%d out of %d printed)\n", NMAX, nz);

	// print error message
	sprintf(m_szerr, "FATAL ERROR: %d zero(s) found on diagonal.", l.size());

}
*/
//=============================================================================
bool NegativeJacobian::m_boutput = false;
bool NegativeJacobian::m_bthrown = false;

//-----------------------------------------------------------------------------
NegativeJacobian::NegativeJacobian(int iel, int ng, double vol, RgElement* pe)
{
	m_iel = iel;
	m_ng = ng;
	m_vol = vol;
	m_pel = pe;
	m_bthrown = true;
	what("Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", m_iel, m_ng + 1, m_vol);
}

//-----------------------------------------------------------------------------
bool NegativeJacobian::DoOutput()
{
	return m_boutput;
}

//-----------------------------------------------------------------------------
void NegativeJacobian::clearFlag()
{
	m_bthrown = false;
}

//-----------------------------------------------------------------------------
bool NegativeJacobian::IsThrown()
{
	return m_bthrown;
}

//-----------------------------------------------------------------------------
NANInResidualDetected::NANInResidualDetected(const FENodalDofInfo& ndi)
{
	what("NAN detected in residual vector at index %d.\nNode id = %d, dof = %d ('%s')", ndi.m_eq, ndi.m_node, ndi.m_dof, ndi.szdof);
}

//-----------------------------------------------------------------------------
NANInSolutionDetected::NANInSolutionDetected(const FENodalDofInfo& ndi)
{
	what("NAN detected in solution vector at index %d.\nNode id = %d, dof = %d ('%s')", ndi.m_eq, ndi.m_node, ndi.m_dof, ndi.szdof);
}


//-----------------------------------------------------------------------------
FEMultiScaleException::FEMultiScaleException(int eid, int gpt)
{
	what("The RVE problem has failed at element %d, gauss point %d.\nAborting macro run.", eid, gpt + 1);
}

//-----------------------------------------------------------------------------
ConcentrationChangeDetected::ConcentrationChangeDetected(const FENodalDofInfo& ndi)
{
    what("Sudden concentration change detected.\n");
}

