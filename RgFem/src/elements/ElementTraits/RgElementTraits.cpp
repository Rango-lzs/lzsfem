#include "RgElementTraits.h"

RgElementTraits::RgElementTraits(int ni, int ne, ElementCategory c, ElementShape s, ElementType t)
{
	m_neln = ne;
	m_nint = ni;
	m_faces = 0;
	m_spec.eclass = c;
	m_spec.eshape = s;
	m_spec.etype  = t;
	m_H.resize(ni, ne);
}

RgGaussPoint RgElementTraits::gaussPoint(int i) const
{
	return gaussPoints[i];
}
