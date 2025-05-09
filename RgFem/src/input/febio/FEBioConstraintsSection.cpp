#include "FEBioConstraintsSection.h"
#include "femcore/FEModel.h"
#include "femcore/FESurfaceConstraint.h"
#include "femcore/FENodeSetConstraint.h"
#include "femcore/FESurfacePairConstraintNL.h"
#include "femcore/FEModelLoad.h"
#include "femcore/FEMesh.h"
#include "femcore/FESurface.h"
#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEInitialCondition.h"

void FEBioConstraintsSection25::Parse(XMLTag &tag)
{
	// make sure there is something to read
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	++tag;
	do
	{
		if (tag == "constraint")
		{
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0)
			{
				// check the name attribute
				const char* szname = tag.AttributeValue("name");
				if (szname == 0) throw XMLReader::InvalidAttributeValue(tag, "name", "(unknown)");

				// make sure this is a leaf
				if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

				// see if we can find this constraint
				FEModel& fem = *GetFEModel();
				int NLC = fem.NonlinearConstraints();
				FENLConstraint* plc = 0;
				for (int i=0; i<NLC; ++i)
				{
					FENLConstraint* pci = fem.NonlinearConstraint(i);
					if (pci->GetName() == szname) { plc = pci; }
				}
				if (plc == 0) throw XMLReader::InvalidAttributeValue(tag, "name", szname);

				// add this boundary condition to the current step
				GetBuilder()->AddComponent(plc);
			}
			else
			{
				FENLConstraint* plc = fecore_new<FENLConstraint>(sztype, GetFEModel());
				if (plc == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

				const char* szname = tag.AttributeValue("name", true);
				if (szname) plc->SetName(szname);

				// get the surface
				// Note that not all constraints define a surface
				FESurfaceConstraint* psc = dynamic_cast<FESurfaceConstraint*>(plc);
				if (psc && psc->GetSurface())
				{
					FESurface* psurf = psc->GetSurface();
					mesh.AddSurface(psurf);
					const char* szsurf = tag.AttributeValue("surface");
					FEFacetSet* pface = mesh.FindFacetSet(szsurf);
					if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);
					if (GetBuilder()->BuildSurface(*psurf, *pface, true) == false) throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);
				}

                // get the nodeset for other constraints
                // Note that not all constraints define a nodeset
                FENodeSetConstraint* pns = dynamic_cast<FENodeSetConstraint*>(plc);
                if (pns && pns->GetNodeSet())
                {
                    FENodeSet* pnset = pns->GetNodeSet();
                    const char* sznset = tag.AttributeValue("node_set");
                    FENodeSet* pset = mesh.FindNodeSet(sznset);
                    if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sznset);
                    pnset->Add(pset->GetNodeList());
                }

                // get the surface pair
                FESurfacePairConstraintNL* pspc = dynamic_cast<FESurfacePairConstraintNL*>(plc);
                if (pspc && pspc->GetPrimarySurface() && pspc->GetSecondarySurface())
                {
                    // get the surface pair
                    const char* szpair = tag.AttributeValue("surface_pair");
                    FESurfacePair* surfacePair =mesh.FindSurfacePair(szpair);
                    if (surfacePair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
                    
                    // build the surfaces
                    if (GetBuilder()->BuildSurface(*pspc->GetSecondarySurface(), *surfacePair->GetSecondarySurface(), pspc->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
                    if (GetBuilder()->BuildSurface(*pspc->GetPrimarySurface(), *surfacePair->GetPrimarySurface(), pspc->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

                    // Make sure we have both surfaces
                    FESurface* pss = pspc->GetPrimarySurface (); if ((pss == 0) || (pss->Elements()==0)) throw XMLReader::MissingAttribute(tag,"Missing constraint primary surface");
                    mesh.AddSurface(pss);
                    FESurface* pms = pspc->GetSecondarySurface(); if ((pms == 0) || (pms->Elements()==0)) throw XMLReader::MissingAttribute(tag,"Missing constraint secondary surface");
                    mesh.AddSurface(pms);
                }
                
				// read the parameter list
				ReadParameterList(tag, plc);

				// add this constraint to the current step
				GetBuilder()->AddNonlinearConstraint(plc);
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//--------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioConstraintsSection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	// count nr of faces
	int faces = 0, N, nf[9];
	XMLTag t(tag); ++t;
	while (!t.isend()) { faces++; ++t; }

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
			if      (tag == "quad4") el.SetType(FE_QUAD4NI);
			else if (tag == "tri3" ) el.SetType(FE_TRI3NI );
			else if (tag == "tri6" ) el.SetType(FE_TRI6NI);
            else if (tag == "quad8" ) el.SetType(FE_QUAD8NI);
            else if (tag == "quad9" ) el.SetType(FE_QUAD9NI);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(feb->m_ntri3);
			else if (tag == "tri6" ) el.SetType(feb->m_ntri6);
			else if (tag == "tri7" ) el.SetType(feb->m_ntri7);
			else if (tag == "tri10") el.SetType(feb->m_ntri10);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else if (tag == "quad9") el.SetType(FE_QUAD9G9);
			else throw XMLReader::InvalidTag(tag);
		}

		N = el.Nodes();

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
			FEElement* pe = m.FindElementFromID(nf[0]);
			if (pe)
			{
				int ne[9];
				int nn = pe->GetFace(nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
				el.m_elem[0] = pe;
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}
