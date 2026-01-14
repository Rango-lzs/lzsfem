#include "elements/RgElement/RgShellElement.h"
#include "basicio/DumpStream.h"
using namespace std;

//=================================================================================================
// RgShellElement
//=================================================================================================

RgShellElement::RgShellElement()
{
	m_elem[0] = m_elem[1] = -1;
}

RgShellElement::RgShellElement(const RgShellElement& el)
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
RgShellElement& RgShellElement::operator = (const RgShellElement& el)
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

void RgShellElement::SetTraits(FEElementTraits* ptraits)
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

void RgShellElement::Serialize(DumpStream &ar)
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
// RgShellElementOld
//=================================================================================================
RgShellElementOld::RgShellElementOld()
{
}

RgShellElementOld::RgShellElementOld(const RgShellElementOld& el) : RgShellElement(el)
{
	m_D0 = el.m_D0;
}

//! assignment operator
RgShellElementOld& RgShellElementOld::operator = (const RgShellElementOld& el)
{
	// copy base class
	RgShellElement::operator=(el);

	// copy this class data
	m_D0 = el.m_D0;

	return (*this);
}

void RgShellElementOld::SetTraits(FEElementTraits* ptraits)
{
	RgShellElement::SetTraits(ptraits);
	m_D0.resize(NodeSize());
}

void RgShellElementOld::Serialize(DumpStream& ar)
{
	RgShellElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_D0;
}

//=================================================================================================
// RgShellElementNew
//=================================================================================================

RgShellElementNew::RgShellElementNew()
{

}

RgShellElementNew::RgShellElementNew(const RgShellElementNew& el) : RgShellElement(el)
{
	// TODO: What about all the EAS parameters?
}

//! assignment operator
RgShellElementNew& RgShellElementNew::operator = (const RgShellElementNew& el)
{
	RgShellElement::operator=(el);

	// TODO: What about all the EAS parameters?

	return (*this);
}

void RgShellElementNew::SetTraits(FEElementTraits* ptraits)
{
	RgShellElement::SetTraits(ptraits);

	// TODO: What about all the EAS parameters?
}

void RgShellElementNew::Serialize(DumpStream &ar)
{
	RgShellElement::Serialize(ar);
	ar & m_fa;
	ar & m_Kaai;
	ar & m_alpha;
	ar & m_alphai;
	ar & m_alphat;
	ar & m_Kua;
	ar & m_Kwa;
	ar & m_E;
}

