#include "elements/RgElement.h"
#include "basicio/DumpStream.h"
#include "materials/FEMaterialPoint.h"
#include <math.h>


double* FEElement::H(int order, int n)
{
	if (order == -1) return m_pTraits->m_H[n];
	else return m_pTraits->m_Hp[order][n];
}

int FEElement::ShapeFunctions(int order) const
{
	return (order == -1 ? NodeSize() : m_pTraits->ShapeFunctions(order));
}

int FEElement::GaussPointSize() const {
    return m_pTraits->m_nint;
}

////-----------------------------------------------------------------------------
//bool FEElement::HasNode(int n) const
//{
//	int l = Nodes();
//	for (int i = 0; i<l; ++i)
//		if (m_node[i] == n) return true;
//	return false;
//}
//
////-----------------------------------------------------------------------------
//// see if this element has the list of nodes n. Return 0 if not, 1 if same order
//// and -1 if opposite order
//int FEElement::HasNodes(int* n, const int ns) const
//{
//    int order = 1;
//    int l = Nodes();
//    if (l < ns) return 0;
//    vector<int> num(ns,-1);
//    for (int j=0; j<ns; ++j) {
//        for (int i = 0; i<l; ++i)
//            if (m_node[i] == n[j]) num[j] = i;
//    }
//    for (int j=0; j<ns; ++j) {
//        if (num[j] == -1)
//            return 0;
//    }
//    if ((num[1] - num[0] < 0) || (num[2] - num[1] < 0)) order = -1;
//
//    return order;
//}
//
////-----------------------------------------------------------------------------
//int FEElement::FindNode(int n) const
//{
//	int l = Nodes();
//	for (int i = 0; i<l; ++i)
//		if (m_node[i] == n) return i;
//	return -1;
//}

//-----------------------------------------------------------------------------
FEElement::FEElement() : m_pTraits(0) 
{ 
	static int n = 1;
	m_id = n++;
	/*m_lm = -1;
	m_val = 0.0;
	m_lid = -1;
	m_part = nullptr;
	m_status = ACTIVE;*/
}


//! get the element ID
int FEElement::getId() const { return m_id; }

//! set the element ID
void FEElement::setId(int n) { m_id = n; }

//! Get the element's material ID
int FEElement::getMatId() const { return m_mat_id; }

//! Set the element's material ID
void FEElement::setMatId(int id) { m_mat_id = id; }

void FEElement::setLocalId(int lid)
{
    m_loc_id = lid;
}

int FEElement::getLocalId() const
{
    return m_loc_id;
}

const std::vector<NodeId>& FEElement::getNodeIds()
{
    return m_node;
}

//-----------------------------------------------------------------------------
void FEElement::SetTraits(FEElementTraits* ptraits)
{
	m_pTraits = ptraits;
	m_node.resize(NodeSize());
	//m_lnode.resize(NodeSize());
	m_state.Create(GaussPointSize());
}

//! serialize
void FEElement::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		int type = getType();
		ar << type;
		/*ar << m_nID << m_lid << m_mat;
		ar << m_node;
		ar << m_lnode;
		ar << m_lm << m_val;
		ar << m_status;*/
	}
	else
	{
		int ntype;
		ar >> ntype; setType(ntype);
        /*ar >> m_nID >> m_lid >> m_mat;
        ar >> m_node;
        ar >> m_lnode;
        ar >> m_lm >> m_val;
        ar >> m_status;*/
	}
}

////! return the nodes of the face
//int FEElement::GetFace(int nface, int* nf) const
//{
//	int nn = -1;
//	const int* en = &(m_node[0]);
//	switch (Shape())
//	{
//	case ET_HEX8:
//		nn = 4;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; break;
//		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; break;
//		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; break;
//		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; break;
//		}
//		break;
//	case ET_PENTA6:
//		switch (nface)
//		{
//		case 0: nn = 4; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; break;
//		case 1: nn = 4; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; break;
//		case 2: nn = 4; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; break;
//		case 3: nn = 3; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[1]; break;
//		case 4: nn = 3; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[5]; break;
//		}
//		break;
//	case ET_PENTA15:
//		switch (nface)
//		{
//		case 0: nn = 8; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; nf[4] = en[6]; nf[5] = en[13]; nf[6] = en[9]; nf[7] = en[12]; break;
//		case 1: nn = 8; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[7]; nf[5] = en[14]; nf[6] = en[10]; nf[7] = en[13]; break;
//		case 2: nn = 8; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; nf[4] = en[12]; nf[5] = en[11]; nf[6] = en[14]; nf[7] = en[8]; break;
//		case 3: nn = 6; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[8]; nf[4] = en[7]; nf[5] = en[6]; break;
//		case 4: nn = 6; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[9]; nf[4] = en[10]; nf[5] = en[11]; break;
//		}
//		break;
//	case ET_PYRA5:
//		switch (nface)
//		{
//		case 0: nn = 3; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; break;
//		case 1: nn = 3; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[4]; break;
//		case 2: nn = 3; nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[4]; break;
//		case 3: nn = 3; nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; break;
//		case 4: nn = 4; nf[0] = en[3]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[0]; break;
//		}
//		break;
//	case ET_TET4:
//	case ET_TET5:
//		nn = 3;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = nf[3] = en[3]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = nf[3] = en[3]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = nf[3] = en[3]; break;
//		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = nf[3] = en[0]; break;
//		}
//		break;
//	case ET_TET10:
//		nn = 6;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; break;
//		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; break;
//		}
//		break;
//	case ET_TET15:
//		nn = 7;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; nf[6] = en[11]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; nf[6] = en[12]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; nf[6] = en[13]; break;
//		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; nf[6] = en[10]; break;
//		}
//		break;
//	case ET_TET20:
//		nn = 10;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[5]; nf[5] = en[12]; nf[6] = en[13]; nf[7] = en[10]; nf[8] = en[11]; nf[9] = en[16]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[14]; nf[6] = en[15]; nf[7] = en[13]; nf[8] = en[14]; nf[9] = en[17]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[9]; nf[4] = en[8]; nf[5] = en[10]; nf[6] = en[11]; nf[7] = en[14]; nf[8] = en[15]; nf[9] = en[18]; break;
//		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[7]; nf[4] = en[6]; nf[5] = en[5]; nf[6] = en[4]; nf[7] = en[10]; nf[8] = en[8]; nf[9] = en[19]; break;
//		}
//		break;
//	case ET_HEX20:
//		nn = 8;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; break;
//		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; break;
//		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[9]; nf[7] = en[8]; break;
//		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; break;
//		}
//		break;
//	case ET_HEX27:
//		nn = 9;
//		switch (nface)
//		{
//		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; nf[8] = en[20]; break;
//		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; nf[8] = en[21]; break;
//		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; nf[8] = en[22]; break;
//		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; nf[8] = en[23]; break;
//		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[9]; nf[7] = en[8]; nf[8] = en[24]; break;
//		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; nf[8] = en[25]; break;
//		}
//		break;
//	case ET_QUAD4:
//		nn = 4;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3];
//		break;
//	case ET_QUAD8:
//		nn = 8;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7];
//		break;
//	case ET_QUAD9:
//		nn = 9;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7]; nf[8] = en[8];
//		break;
//	case ET_TRI3:
//		nn = 3;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2];
//		break;
//	case ET_TRI6:
//		nn = 6;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5];
//		break;
//	case ET_TRI7:
//		nn = 7;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6];
//		break;
//	case ET_TRI10:
//		nn = 7;
//		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7]; nf[8] = en[8]; nf[9] = en[9];
//		break;
//	}
//
//	return nn;
//}


double* FEElement::H(int n)
{
    return m_pTraits->m_H[n];
}

const double* FEElement::H(int n) const
{
    return m_pTraits->m_H[n];
}

//! Get the material point data
FEMaterialPoint* FEElement::GetMaterialPoint(int n)
{
    return m_state[n];
}

//! set the material point data
void FEElement::SetMaterialPointData(FEMaterialPoint* pmp, int n)
{
      pmp->m_elem = this;
      pmp->m_index = n;
      m_state[n] = pmp;
 }




////-----------------------------------------------------------------------------
//FEDiscreteElement::FEDiscreteElement(const FEDiscreteElement& el)
//{
//	// set the traits of the element
//	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }
//
//	// copy base class data
//	m_mat = el.m_mat;
//	m_nID = el.m_nID;
//	m_lid = el.m_lid;
//	m_node = el.m_node;
//	m_lnode = el.m_lnode;
//	m_lm = el.m_lm;
//	m_val = el.m_val;
//}
//
//FEDiscreteElement& FEDiscreteElement::operator =(const FEDiscreteElement& el)
//{
//	// set the traits of the element
//	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }
//
//	// copy base class data
//	m_mat = el.m_mat;
//	m_nID = el.m_nID;
//	m_lid = el.m_lid;
//	m_node = el.m_node;
//	m_lnode = el.m_lnode;
//	m_lm = el.m_lm;
//	m_val = el.m_val;
//
//	return (*this);
//}
//
////-----------------------------------------------------------------------------
//FELineElement::FELineElement()
//{
//	m_lid = -1;
//}
//
//FELineElement::FELineElement(const FELineElement& el)
//{
//	// set the traits of the element
//	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }
//
//	// copy data
//	m_lid = el.m_lid;
//
//	// copy base class data
//	m_mat = el.m_mat;
//	m_nID = el.m_nID;
//	m_lid = el.m_lid;
//	m_node = el.m_node;
//	m_lnode = el.m_lnode;
//	m_lm = el.m_lm;
//	m_val = el.m_val;
//}
//
//FELineElement& FELineElement::operator = (const FELineElement& el)
//{
//	// set the traits of the element
//	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }
//
//	// copy data
//	m_lid = el.m_lid;
//
//	// copy base class data
//	m_mat = el.m_mat;
//	m_nID = el.m_nID;
//	m_lid = el.m_lid;
//	m_node = el.m_node;
//	m_lnode = el.m_lnode;
//	m_lm = el.m_lm;
//	m_val = el.m_val;
//
//	return (*this);
//}
//
//void FELineElement::SetTraits(FEElementTraits* pt)
//{
//	// we don't allocate state data for surface elements
//	m_pTraits = pt;
//	m_node.resize(Nodes());
//	m_lnode.resize(Nodes());
//}
