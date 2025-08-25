#include "elements/FEShellElement.h"
#include "basicio/DumpStream.h"
using namespace std;

//=================================================================================================
// FEShellElement
//=================================================================================================

FEShellElement::FEShellElement()
{
	m_elem[0] = m_elem[1] = -1;
}

FEShellElement::FEShellElement(const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }

	// copy base class data

	// copy shell data
	m_h0 = el.m_h0;
	m_ht = el.m_ht;
	m_d0 = el.m_d0;
	m_g0[0] = el.m_g0[0]; m_g0[1] = el.m_g0[1]; m_g0[2] = el.m_g0[2];
	m_gt[0] = el.m_gt[0]; m_gt[1] = el.m_gt[1]; m_gt[2] = el.m_gt[2];
	m_gp[0] = el.m_gp[0]; m_gp[1] = el.m_gp[1]; m_gp[2] = el.m_gp[2];

	m_G0[0] = el.m_G0[0]; m_G0[1] = el.m_G0[1]; m_G0[2] = el.m_G0[2];
	m_Gt[0] = el.m_Gt[0]; m_Gt[1] = el.m_Gt[1]; m_Gt[2] = el.m_Gt[2];

	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];
}

//! assignment operator
FEShellElement& FEShellElement::operator = (const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pTraits) { SetTraits(el.m_pTraits); m_state = el.m_state; }

	// copy base class data
	
	// copy shell data
	m_h0 = el.m_h0;
	m_ht = el.m_ht;
	m_d0 = el.m_d0;
	m_g0[0] = el.m_g0[0]; m_g0[1] = el.m_g0[1]; m_g0[2] = el.m_g0[2];
	m_gt[0] = el.m_gt[0]; m_gt[1] = el.m_gt[1]; m_gt[2] = el.m_gt[2];
	m_gp[0] = el.m_gp[0]; m_gp[1] = el.m_gp[1]; m_gp[2] = el.m_gp[2];

	m_G0[0] = el.m_G0[0]; m_G0[1] = el.m_G0[1]; m_G0[2] = el.m_G0[2];
	m_Gt[0] = el.m_Gt[0]; m_Gt[1] = el.m_Gt[1]; m_Gt[2] = el.m_Gt[2];

	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

void FEShellElement::SetTraits(FEElementTraits* ptraits)
{
	FEElement::SetTraits(ptraits);
	m_h0.assign(NodeSize(), 0.0);
	m_ht.assign(NodeSize(), 0.0);
	m_d0.assign(NodeSize(), Vector3d(0, 0, 0));
	m_g0[0].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_g0[1].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_g0[2].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gt[0].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gt[1].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gt[2].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gp[0].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gp[1].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_gp[2].assign(GaussPointSize(), Vector3d(0, 0, 0));

	m_G0[0].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_G0[1].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_G0[2].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_Gt[0].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_Gt[1].assign(GaussPointSize(), Vector3d(0, 0, 0));
	m_Gt[2].assign(GaussPointSize(), Vector3d(0, 0, 0));
}

void FEShellElement::Serialize(DumpStream &ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false) {
		ar & m_h0;
		ar & m_d0;
		ar & m_g0[0] & m_g0[1] & m_g0[2];
		ar & m_gt[0] & m_gt[1] & m_gt[2];
		ar & m_gp[0] & m_gp[1] & m_gp[2];
		ar & m_G0[0] & m_G0[1] & m_G0[2];
		ar & m_Gt[0] & m_Gt[1] & m_Gt[2];
		ar & m_elem[0] & m_elem[1];
	}
	ar & m_ht;
}

//=================================================================================================
// FEShellElementOld
//=================================================================================================
FEShellElementOld::FEShellElementOld()
{
}

FEShellElementOld::FEShellElementOld(const FEShellElementOld& el) : FEShellElement(el)
{
	m_D0 = el.m_D0;
}

//! assignment operator
FEShellElementOld& FEShellElementOld::operator = (const FEShellElementOld& el)
{
	// copy base class
	FEShellElement::operator=(el);

	// copy this class data
	m_D0 = el.m_D0;

	return (*this);
}

void FEShellElementOld::SetTraits(FEElementTraits* ptraits)
{
	FEShellElement::SetTraits(ptraits);
	m_D0.resize(NodeSize());
}

void FEShellElementOld::Serialize(DumpStream& ar)
{
	FEShellElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_D0;
}

//=================================================================================================
// FEShellElementNew
//=================================================================================================

FEShellElementNew::FEShellElementNew()
{

}

FEShellElementNew::FEShellElementNew(const FEShellElementNew& el) : FEShellElement(el)
{
	// TODO: What about all the EAS parameters?
}

//! assignment operator
FEShellElementNew& FEShellElementNew::operator = (const FEShellElementNew& el)
{
	FEShellElement::operator=(el);

	// TODO: What about all the EAS parameters?

	return (*this);
}

void FEShellElementNew::SetTraits(FEElementTraits* ptraits)
{
	FEShellElement::SetTraits(ptraits);

	// TODO: What about all the EAS parameters?
}

void FEShellElementNew::Serialize(DumpStream &ar)
{
	FEShellElement::Serialize(ar);
	ar & m_fa;
	ar & m_Kaai;
	ar & m_alpha;
	ar & m_alphai;
	ar & m_alphat;
	ar & m_Kua;
	ar & m_Kwa;
	ar & m_E;
}
