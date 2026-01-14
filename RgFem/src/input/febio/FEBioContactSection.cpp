#include "FEBioContactSection.h"
//#include "femcore/FEAugLagLinearConstraint.h"
#include "femcore/FENodalBC.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"

//-----------------------------------------------------------------------------
FEBioContactSection::MissingPrimarySurface::MissingPrimarySurface()
{
	SetErrorString("Missing contact primary surface");
}

//-----------------------------------------------------------------------------
FEBioContactSection::MissingSecondarySurface::MissingSecondarySurface()
{
	SetErrorString("Missing contact secondary surface");
}

//-----------------------------------------------------------------------------
// --- L I N E A R   C O N S T R A I N T ---
void FEBioContactSection::ParseLinearConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
    DOFS& dofs = fem.GetDOFS();
	FEMesh& m = fem.GetMesh();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	//// create a new linear constraint manager
	//FELinearConstraintSet* pLCS = dynamic_cast<FELinearConstraintSet*>(fecore_new<FENLConstraint>("linear constraint", GetFEModel()));
	//fem.AddNonlinearConstraint(pLCS);

	// read the linear constraints
	++tag;
	do
	{
		//if (tag == "linear_constraint")
		//{
		//	FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, &fem);

		//	++tag;
		//	do
		//	{
		//		int node, bc;
		//		double val;
		//		if (tag == "node")
		//		{
		//			tag.value(val);
		//			
		//			const char* szid = tag.AttributeValue("id");
		//			node = atoi(szid);

		//			const char* szbc = tag.AttributeValue("bc");
  //                  int ndof = dofs.GetDOF(szbc);
  //                  if (ndof >= 0) bc = ndof;
		//			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

		//			pLC->AddDOF(node, bc, val);
		//		}
		//		else throw XMLReader::InvalidTag(tag);
		//		++tag;
		//	}
		//	while (!tag.isend());

		//	// add the linear constraint to the system
		//	pLCS->add(pLC);
		//}
		//else if (ReadParameter(tag, pLCS) == false)
		//{
		//	throw XMLReader::InvalidTag(tag);
		//}
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioContactSection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	int N, nf[9];

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	s.Create(faces);

	FEModelBuilder* feb = GetBuilder();

	// read faces
	++tag;
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		// set the element type/integration rule
		if (bnodal)
		{
			if      (tag == "quad4") el.setType(FE_QUAD4NI);
			else if (tag == "tri3" ) el.setType(FE_TRI3NI );
			else if (tag == "tri6" ) el.setType(FE_TRI6NI );
            else if (tag == "quad8" ) el.setType(FE_QUAD8NI);
            else if (tag == "quad9" ) el.setType(FE_QUAD9NI);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.setType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.setType(feb->m_ntri3);
			else if (tag == "tri6" ) el.setType(feb->m_ntri6);
			else if (tag == "tri7" ) el.setType(feb->m_ntri7);
			else if (tag == "tri10") el.setType(feb->m_ntri10);
			else if (tag == "quad8") el.setType(FE_QUAD8G9);
			else if (tag == "quad9") el.setType(FE_QUAD9G9);
			else throw XMLReader::InvalidTag(tag);
		}

		N = el.NodeSize();

		if (nfmt == 0)
		{
			tag.value(nf, N);
			for (int j=0; j<N; ++j) 
			{
				int nid = nf[j]-1;
				if ((nid<0)||(nid>= NN)) throw XMLReader::InvalidValue(tag);
				el.m_node[j] = nid;
			}
		}
		else if (nfmt == 1)
		{
			tag.value(nf, 2);
			RgElement* pe = m.FindElementFromID(nf[0]);
			if (pe)
			{
				int ne[4];
				int nn = pe->GetFace(nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
				el.m_elem[0] = pe;
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}

	s.InitSurface();
	s.CreateRgMaterialPointData();

	return true;
}
//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection4::Parse(XMLTag& tag)
{
	// make sure there are children
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();

	// loop over tags
	++tag;
	do
	{
		if (tag == "contact")
		{
			// get the contact type
			const char* sztype = tag.AttributeValue("type");

			//// Try to create a contact interface using the FEBio kernel. 
			//FESurfacePairConstraint* pci = fecore_new<FESurfacePairConstraint>(sztype, &fem);
			//if (pci == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			/*GetBuilder()->AddContactInterface(pci);
			ParseContactInterface(tag, pci);*/
		}

		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioContactSection4::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	//// get the surface pair
	//const char* szpair = tag.AttributeValue("surface_pair");
	//FESurfacePair* surfacePair = m.FindSurfacePair(szpair);
	//if (surfacePair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	//// build the surfaces
	//if (GetBuilder()->BuildSurface(*pci->GetSecondarySurface(), *surfacePair->GetSecondarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
	//if (GetBuilder()->BuildSurface(*pci->GetPrimarySurface(), *surfacePair->GetPrimarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	//// get the parameter list
	//ReadParameterList(tag, pci);

	//// Make sure we have both surfaces
	//FESurface* pss = pci->GetPrimarySurface(); if ((pss == 0) || (pss->Elements() == 0)) throw MissingPrimarySurface();
	//m.AddSurface(pss);
	//FESurface* pms = pci->GetSecondarySurface(); if ((pms == 0) || (pms->Elements() == 0)) throw MissingSecondarySurface();
	//m.AddSurface(pms);
}
