#include "FEBioBoundarySection.h"
#include"femcore/FEModel.h"
//#include"femcore/Domain/FEDiscreteDomain.h"
//#include"femcore/FEAugLagLinearConstraint.h"
#include"femcore/FEPrescribedDOF.h"
#include"femcore/FEFixedBC.h"
#include"femcore/FELinearConstraintManager.h"
//#include"femcore/FEPeriodicLinearConstraint.h"
//#include <FEBioRVE/FEPeriodicLinearConstraint2O.h"
//#include"femcore/FEMergedConstraint.h"
#include"femcore/FEFacetSet.h"
#include"logger/log.h"
#include"femcore/FEModelLoad.h"
#include"femcore/FEInitialCondition.h"

//---------------------------------------------------------------------------------
void FEBioBoundarySection::BuildNodeSetMap()
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	m_NodeSet.clear();
	for (int i = 0; i<m.NodeSets(); ++i)
	{
		FENodeSet* nsi = m.NodeSet(i);
		m_NodeSet[nsi->GetName()] = nsi;
	}
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioBoundarySection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
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
				el.setNodeId(j, nid);
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
				for (int j=0; j<N; ++j) el.setNodeId(j, ne[j]);
				el.m_elem[0] = pe;
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::AddFixedBC(FENodeSet* set, int dof)
{
	FEModel* fem = GetFEModel();
	FEFixedBC* bc = new FEFixedBC(fem);
	bc->SetDOFList(dof);
	bc->SetNodeSet(set);
	fem->AddBoundaryCondition(bc);
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCFix(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();

	// get the DOF indices
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_U = fem.GetDOFIndex("u");
	const int dof_V = fem.GetDOFIndex("v");
	const int dof_W = fem.GetDOFIndex("w");
    const int dof_SX = fem.GetDOFIndex("sx");
    const int dof_SY = fem.GetDOFIndex("sy");
    const int dof_SZ = fem.GetDOFIndex("sz");
    const int dof_WX = fem.GetDOFIndex("wx");
    const int dof_WY = fem.GetDOFIndex("wy");
    const int dof_WZ = fem.GetDOFIndex("wz");
    const int dof_GX = fem.GetDOFIndex("gx");
    const int dof_GY = fem.GetDOFIndex("gy");
    const int dof_GZ = fem.GetDOFIndex("gz");

	// see if a set is defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// read the set
		FEMesh& mesh = fem.GetMesh();
		FENodeSet* ps = mesh.FindNodeSet(szset);
		if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// get the bc attribute
		const char* sz = tag.AttributeValue("bc");

		// Make sure this is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// loop over all nodes in the nodeset
		int ndof = dofs.GetDOF(sz);
		if (ndof >= 0) AddFixedBC(ps, ndof);
		else
		{
			// The supported fixed BC strings don't quite follow the dof naming convention.
			// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
			if      (strcmp(sz, "xy"  ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Y); }
			else if (strcmp(sz, "yz"  ) == 0) { AddFixedBC(ps, dof_Y ); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "xz"  ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "xyz" ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Y); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "uv"  ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_V); }
			else if (strcmp(sz, "vw"  ) == 0) { AddFixedBC(ps, dof_V ); AddFixedBC(ps, dof_W); }
			else if (strcmp(sz, "uw"  ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_W); }
			else if (strcmp(sz, "uvw" ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_V); AddFixedBC(ps, dof_W); }
            else if (strcmp(sz, "sxy" ) == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SY); }
            else if (strcmp(sz, "syz" ) == 0) { AddFixedBC(ps, dof_SY); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "sxz" ) == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "sxyz") == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SY); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "wxy" ) == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WY); }
            else if (strcmp(sz, "wyz" ) == 0) { AddFixedBC(ps, dof_WY); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "wxz" ) == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "wxyz") == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WY); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "gxy" ) == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GY); }
            else if (strcmp(sz, "gyz" ) == 0) { AddFixedBC(ps, dof_GY); AddFixedBC(ps, dof_GZ); }
            else if (strcmp(sz, "gxz" ) == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GZ); }
            else if (strcmp(sz, "gxyz") == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GY); AddFixedBC(ps, dof_GZ); }
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
		}
	}
	else
	{
		// The format where the bc can be defined on each line is no longer supported.
		throw XMLReader::MissingAttribute(tag, "set");
	}
}


//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCPrescribe(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	// see if this tag defines a set
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// Find the set
		FENodeSet* ps = mesh.FindNodeSet(szset);
		if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// get the bc attribute
		const char* sz = tag.AttributeValue("bc");

		// find the dof index from its symbol
		int bc = dofs.GetDOF(sz);
		if (bc == -1) 
		{
			// the temperature degree of freedom was renamed
			// for backward compatibility we need to check for it
			if (strcmp(sz, "t") == 0) bc = dofs.GetDOF("T");
			throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
		}

		// get the lc attribute
		sz = tag.AttributeValue("lc");
		int lc = atoi(sz);

		// make sure this tag is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// get the scale factor
		double s = 1;
		value(tag, s);

		// create the bc
        FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(RANGO_NEW<FEBoundaryCondition>(&fem ,"prescribe"));
		pdc->SetScale(s).SetDOF(bc);

		if (lc >= 0)
		{
			FEParam* p = pdc->GetParameter("scale");
			if (p == nullptr) throw XMLReader::InvalidTag(tag);
			fem.AttachLoadController(p, lc);
		}

		// add this boundary condition to the current step
		GetBuilder()->AddBC(pdc);

		// add nodes in the nodeset
		pdc->SetNodeSet(ps);
	}
	else
	{
		// count how many prescibed nodes there are
		int ndis = tag.children();

		// determine whether prescribed BC is relative or absolute
		bool br = false;
		const char* sztype = tag.AttributeValue("type",true);
		if (sztype && strcmp(sztype, "relative") == 0) br = true;

		// read the prescribed data
		++tag;
		for (int i=0; i<ndis; ++i)
		{
			int n = ReadNodeID(tag);
			const char* sz = tag.AttributeValue("bc");

			// get the dof index from its symbol
			int bc = dofs.GetDOF(sz);
			if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			sz = tag.AttributeValue("lc");
			int lc = atoi(sz)-1;

			double scale;
			tag.value(scale);

			FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(RANGO_NEW<FEBoundaryCondition>(&fem ,"prescribe"));
			pdc->SetDOF(bc);
			pdc->SetScale(scale);
			pdc->SetRelativeFlag(br);

			FENodeSet* ps = new FENodeSet(&fem);
			ps->Add(n);
			pdc->SetNodeSet(ps);

			if (lc >= 0)
			{
				FEParam* p = pdc->GetParameter("scale");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				fem.AttachLoadController(p, lc);
			}

			// add this boundary condition to the current step
			GetBuilder()->AddBC(pdc);

			// next tag
			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseSpringSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	//Rango TODO:
	//// determine the spring type
	//const char* szt = tag.AttributeValue("type", true);
	//if (szt == 0) szt = "linear";
	//FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(szt, &fem));
	//if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	//// create a new spring "domain"
	//FECoreKernel& febio = FECoreKernel::GetInstance();
	//FE_Element_Spec spec;
	//spec.eclass = FE_ELEM_TRUSS;
	//spec.eshape = ET_TRUSS2;
	//spec.etype  = FE_DISCRETE;
	//FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(spec, &mesh, pm));
	//mesh.AddDomain(pd);

	//pd->Create(1, spec);
	//FEDiscreteElement& de = pd->Element(0);
	//de.SetID(++GetBuilder()->m_maxid);
	//
	//// add a new material for each spring
	//fem.AddMaterial(pm);
	//pm->SetID(fem.Materials());
	//pd->SetMatID(fem.Materials()-1);

	// read spring discrete elements
	//++tag;
	//do
	//{
	//	// read the required node tag
	//	if (tag == "node")
	//	{
	//		int n[2];
	//		tag.value(n, 2);
	//		de.m_node[0] = n[0]-1;
	//		de.m_node[1] = n[1]-1;
	//	}
	//	else
	//	{
	//		// read the actual spring material parameters
	//		FEParameterList& pl = pm->GetParameterList();
	//		if (ReadParameter(tag, pl) == 0)
	//		{
	//			throw XMLReader::InvalidTag(tag);
	//		}
	//	}
	//	++tag;
	//}
	//while (!tag.isend());

	//pd->CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
//! Parse the linear constraints section of the xml input file
//! This section is a subsection of the Boundary section

void FEBioBoundarySection::ParseConstraints(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	FEModelBuilder* feb = GetBuilder();

	// read the parent node
	int nodeID;
	tag.AttributeValue("node", nodeID);
	int parentNode = feb->FindNodeFromID(nodeID);

	// get the dofs
	const char* szbc = tag.AttributeValue("bc");
	std::vector<int> dofList;
	dofs.ParseDOFString(szbc, dofList);
	int ndofs = (int) dofList.size();

	// allocate linear constraints
	std::vector<FELinearConstraint*> LC;
	for (int i=0; i<ndofs; ++i)
	{
		int dof = dofList[i];
		if (dof < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

		LC[i] = RANGO_NEW<FELinearConstraint>(&fem, "");
		LC[i]->SetParentDof(dof, parentNode);
	}

	// read the child nodes
	++tag;
	do
	{
		if (tag == "node")
		{
			// get the node
			int childNode = ReadNodeID(tag);

			// get the dof
			// (if ommitted we take the parent dof)
			int childDOF = -1;
			const char* szbc = tag.AttributeValue("bc", true);
			if (szbc)
			{
				childDOF = dofs.GetDOF(szbc);
				if (childDOF < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
			}

			// get the coefficient
			double val;
			tag.value(val);

			// add it to the list
			for (int i=0; i<ndofs; ++i)
			{
				int ndof = (childDOF < 0 ? LC[i]->GetParentDof() : childDOF);
				LC[i]->AddChildDof(ndof, childNode, val);
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add the linear constraint to the system
	for (int i=0; i<ndofs; ++i)
	{
		fem.GetLinearConstraintManager().AddLinearConstraint(LC[i]);
		GetBuilder()->AddComponent(LC[i]);
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();

	// read the parameters
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype = 0;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FESurface& s = *(ntype == 1? pci->GetSecondarySurface() : pci->GetPrimarySurface());
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, pci->UseNodalIntegration());
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}


//-----------------------------------------------------------------------------
//! Parses the contact section of the xml input file
//! The contact section is a subsection of the boundary section

void FEBioBoundarySection::ParseContactSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FEModelBuilder* feb = GetBuilder();

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	// Not all contact interfaces can be parsed automatically.
	// First, check all these special cases.
	if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FESurfacePairConstraint* ps = RANGO_NEW<FESurfacePairConstraint>(GetFEModel() ,szt);
		if (ps)
		{
			fem.AddSurfacePairConstraint(ps);

			++tag;
			do
			{
				if (ReadParameter(tag, ps) == false)
				{
					if (tag == "surface")
					{
						FESurface& s = *ps->GetPrimarySurface();

						int nfmt = 0;
						const char* szfmt = tag.AttributeValue("format", true);
						if (szfmt)
						{
							if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
							else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
						}

						// read the surface section
						ParseSurfaceSection(tag, s, nfmt, true);
					}
					else throw XMLReader::InvalidTag(tag);
				}
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}
	else if (strcmp(szt, "rigid") == 0)
	{
		// --- R I G I D   B O D Y   I N T E R F A C E ---

		// count how many rigid nodes there are
		int nrn= 0;
		XMLTag t(tag); ++t;
		while (!t.isend()) { nrn++; ++t; }

		++tag;
		int id, rb, rbp = -1;
		FENodalBC* prn = 0;
		FENodeSet* ns = 0;
		for (int i=0; i<nrn; ++i)
		{
			id = atoi(tag.AttributeValue("id"))-1;
			rb = atoi(tag.AttributeValue("rb"));

			if ((prn == 0) || (rb != rbp))
			{
                prn = RANGO_NEW<FENodalBC>(&fem ,"FERigidNodeSet");

				prn->SetParameter("rb", rb);

				ns = new FENodeSet(&fem);
				prn->SetNodeSet(ns);

				// the default shell bc depends on the shell formulation
				// hinged shell = 0
				// clamped shell = 1
				prn->SetParameter("clamp_shells", feb->m_default_shell == OLD_SHELL ? 0 : 1);

				feb->AddBC(prn);
				rbp = rb;
			}
			ns->Add(id);

			++tag;
		}
	}
	else if (strcmp(szt, "linear constraint") == 0)
	{
		FEModel& fem = *GetFEModel();

		// make sure there is a constraint defined
		if (tag.isleaf()) return;

		// create a new linear constraint manager
        FELinearConstraintSet* pLCS =
            dynamic_cast<FELinearConstraintSet*>(RANGO_NEW<FENLConstraint>(GetFEModel() ,szt ));
		fem.AddNonlinearConstraint(pLCS);

		// read the linear constraints
		++tag;
		do
		{
			if (tag == "linear_constraint")
			{
                FEAugLagLinearConstraint* pLC = RANGO_NEW<FEAugLagLinearConstraint>(&fem, "");

				++tag;
				do
				{
					int node, bc;
					double val;
					if (tag == "node")
					{
						tag.value(val);

						const char* szid = tag.AttributeValue("id");
						node = atoi(szid);

						const char* szbc = tag.AttributeValue("bc");
						if      (strcmp(szbc, "x") == 0) bc = 0;
						else if (strcmp(szbc, "y") == 0) bc = 1;
						else if (strcmp(szbc, "z") == 0) bc = 2;
						else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

						pLC->AddDOF(node, bc, val);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());

				// add the linear constraint to the system
				pLCS->add(pLC);
			}
			else if (ReadParameter(tag, pLCS) == false)
			{
				throw XMLReader::InvalidTag(tag);
			}
			++tag;
		}
		while (!tag.isend());
	}
	else
	{
		// If we get here, we try to create a contact interface
		// using the FEBio kernel. 
		FESurfacePairConstraint* pci = RANGO_NEW<FESurfacePairConstraint>(GetFEModel() , szt);
		if (pci)
		{
			// add it to the model
			GetBuilder()->AddContactInterface(pci);

			// parse the interface
			ParseContactInterface(tag, pci);
		}
		else 
		{
			// some nonlinear constraints are also defined in the Contact section, so let's try that next.
			// TODO: These are mostly rigid constraints. Therefore, I would like to move this elsewhere (maybe in the new Rigid section?)
            FENLConstraint* pnlc = RANGO_NEW<FENLConstraint>(GetFEModel() ,szt);
			if (pnlc)
			{
				ReadParameterList(tag, pnlc);
				fem.AddNonlinearConstraint(pnlc);
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
	}
}
