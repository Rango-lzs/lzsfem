#pragma once
#include "femcore/FEParam.h"
#include "RgDomain.h"

//-----------------------------------------------------------------------------
class RgDomainParameter : public FEParam
{
public:
	RgDomainParameter(RgDomain* pdom = 0);

	RgDomain* GetDomain() { return m_pdom; }

	void SetDomain(RgDomain* pdom) { m_pdom = pdom; }

private:
	RgDomain*	m_pdom;
};