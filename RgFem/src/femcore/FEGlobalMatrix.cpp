#include "FEGlobalMatrix.h"
#include "FEModel.h"
#include "femcore/Domain/FEDomain.h"
#include "FESurface.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEElementMatrix::FEElementMatrix(const FEElement& el)
{
	m_node = el.m_node;
}

//-----------------------------------------------------------------------------
FEElementMatrix::FEElementMatrix(const FEElementMatrix& ke) : Matrix(ke)
{
	m_node = ke.m_node;
	m_lmi = ke.m_lmi;
	m_lmj = ke.m_lmj;
}

//-----------------------------------------------------------------------------
FEElementMatrix::FEElementMatrix(const FEElementMatrix& ke, double scale)
{
	m_node = ke.m_node;
	m_lmi = ke.m_lmi;
	m_lmj = ke.m_lmj;
	Matrix& T = *this;
	const Matrix& K = ke;
	T = (scale == 1.0 ? K : K*scale);
}

//-----------------------------------------------------------------------------
FEElementMatrix::FEElementMatrix(const FEElement& el, const std::vector<int>& lmi) : Matrix((int)lmi.size(), (int)lmi.size())
{
	m_node = el.m_node;
	m_lmi = lmi;
	m_lmj = lmi;
}

//-----------------------------------------------------------------------------
FEElementMatrix::FEElementMatrix(const FEElement& el, std::vector<int>& lmi, std::vector<int>& lmj) : Matrix((int)lmi.size(), (int)lmj.size())
{
	m_node = el.m_node;
	m_lmi = lmi;
	m_lmj = lmj;
};

//-----------------------------------------------------------------------------
// assignment operator
void FEElementMatrix::operator = (const Matrix& ke)
{
	Matrix::operator=(ke);
}


//-----------------------------------------------------------------------------
//! Takes a SparseMatrix structure that defines the structure of the global Matrix.
FEGlobalMatrix::FEGlobalMatrix(SparseMatrix* pK, bool del)
{
	m_pA = pK;
	m_LM.resize(MAX_LM_SIZE);
	m_pMP = 0;
	m_nlm = 0;
	m_delA = del;
}

//-----------------------------------------------------------------------------
//! Deletes the SparseMatrix variable.
FEGlobalMatrix::~FEGlobalMatrix()
{
	if (m_delA) delete m_pA;
	m_pA = 0;
	if (m_pMP) delete m_pMP;
}

//-----------------------------------------------------------------------------
void FEGlobalMatrix::Clear()
{ 
	if (m_pA) m_pA->Clear(); 
}

//-----------------------------------------------------------------------------
//! Start building the profile. That is delete the old profile (if there was one)
//! and create a new one. 
void FEGlobalMatrix::build_begin(int neq)
{
	if (m_pMP) delete m_pMP;
	m_pMP = new SparseMatrixProfile(neq, neq);

	// initialize it to a diagonal Matrix
	// TODO: Is this necessary?
	m_pMP->CreateDiagonal();

	m_nlm = 0;
}

//-----------------------------------------------------------------------------
//! Add an "element" to the Matrix profile. The definition of an element is quite
//! general at this point. Any two or more degrees of freedom that are connected 
//! somehow are considered an element. This function takes one argument, namely a
//! list of degrees of freedom that are part of such an "element". Elements are 
//! first copied into a local LM array.  When this array is full it is flushed and
//! all elements in this array are added to the Matrix profile.
void FEGlobalMatrix::build_add(std::vector<int>& lm)
{
	if (lm.empty() == false)
	{
		m_LM[m_nlm++] = lm;
		if (m_nlm >= MAX_LM_SIZE) build_flush();
	}
}

//-----------------------------------------------------------------------------
//! Flush the LM array. The LM array stores a buffer of elements that have to be
//! added to the profile. When this buffer is full it needs to be flushed. This
//! flushin operation causes the actual update of the Matrix profile.
void FEGlobalMatrix::build_flush()
{
	int i, j, n, *lm;

	// Since prescribed dofs have an equation number of < -1 we need to modify that
	// otherwise no storage will be allocated for these dofs (even not diagonal elements!).
	for (i=0; i<m_nlm; ++i)
	{
		n = (int)m_LM[i].size();
		if (n > 0)
		{
			lm = &(m_LM[i])[0];
			for (j=0; j<n; ++j) if (lm[j] < -1) lm[j] = -lm[j]-2;
		}
	}

	m_pMP->UpdateProfile(m_LM, m_nlm);
	m_nlm = 0;
}

//-----------------------------------------------------------------------------
//! This function makes sure the LM buffer is flushed and creates the actual
//! sparse Matrix from the Matrix profile.
void FEGlobalMatrix::build_end()
{
	if (m_nlm > 0) build_flush();
	m_pA->Create(*m_pMP);
}

//-----------------------------------------------------------------------------
bool FEGlobalMatrix::Create(FEModel* pfem, int neq, bool breset)
{
	// The first time we come here we build the "static" profile.
	// This static profile stores the contribution to the Matrix profile
	// of the "elements" that do not change. Most elements are static except
	// for instance contact elements which can change connectivity in between
	// calls to the Create() function. Storing the static profile instead of
	// reconstructing it every time we come here saves us a lot of time. The 
	// static profile is stored in the variable m_MPs.

	// begin building the profile
	build_begin(neq);
	{
		// The first time we are here we construct the "static"
		// profile. This profile contains the contribution from
		// all static elements. A static element is defined as
		// an element that never changes its connectity. This 
		// static profile is stored in the MP object. Next time
		// we come here we simply copy the MP object in stead
		// of building it from scratch.
		if (breset)
		{
			m_MPs.Clear();

			// build the Matrix profile
			pfem->BuildMatrixProfile(*this, true);

			// copy the static profile to the MP object
			// Make sure the LM buffer is flushed first.
			build_flush();
			m_MPs = *m_pMP;
		}
		else
		{
			// copy the old static profile
			*m_pMP = m_MPs;
		}

		// Add the "dynamic" profile
		pfem->BuildMatrixProfile(*this, false);
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse Matrix. This is done in the following function
	build_end();

	return true;
}

//-----------------------------------------------------------------------------
//! Constructs the stiffness Matrix from a FEMesh object. 
bool FEGlobalMatrix::Create(FEMesh& mesh, int neq)
{
	// begin building the profile
	build_begin(neq);
	{
		// Add all elements to the profile
		// Loop over all active domains
		for (int nd=0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);
			d.BuildMatrixProfile(*this);
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse Matrix. This is done in the following function
	build_end();

	return true;
}

//-----------------------------------------------------------------------------
//! construct the stiffness Matrix from a mesh
bool FEGlobalMatrix::Create(FEMesh& mesh, int nstart, int nend)
{
	if (nstart > nend) return false;
	int neq = nend - nstart + 1;

	// begin building the profile
	build_begin(neq);
	{
		// Add all elements to the profile
		// Loop over all active domains
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);

			std::vector<int> elm, lm;
			const int NE = d.Elements();
			for (int j = 0; j<NE; ++j)
			{
				FEElement& el = d.ElementRef(j);
				d.UnpackLM(el, elm);

				lm.clear();
				for (int k = 0; k < elm.size(); ++k)
				{
					int nk = elm[k];
					if (nk < -1) nk = -nk - 2;

					if ((nk >= nstart) && (nk <= nend))
					{
						if (elm[k] < -1) lm.push_back(elm[k] + nstart);
						else lm.push_back(elm[k] - nstart);
					}
				}

				build_add(lm);
			}
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse Matrix. This is done in the following function
	build_end();

	return true;
}

//! construct a stiffness Matrix from a surface
bool FEGlobalMatrix::Create(const FESurface& surf, const std::vector<int>& equationIDs)
{
	int N = surf.Nodes();
	if ((int)equationIDs.size() != N) return false;

	// count equations
	int neq = 0;
	for (int i = 0; i < N; ++i) if (equationIDs[i] != -1) neq++;
	if (neq == 0) return false;

	// build the Matrix, assuming one degree of freedom per node
	build_begin(neq);
	for (int i = 0; i<surf.Elements(); ++i) {
		const FESurfaceElement& el = surf.Element(i);
		std::vector<int> elm(el.NodeSize(), -1);
		for (int j = 0; j<el.NodeSize(); ++j)
			elm[j] = equationIDs[el.m_loc_node[j]];
		build_add(elm);
	}
	build_end();

	return true;
}

void FEGlobalMatrix::Assemble(const FEElementMatrix& ke)
{
	m_pA->Assemble(ke, ke.RowIndices(), ke.ColumnsIndices());
}
