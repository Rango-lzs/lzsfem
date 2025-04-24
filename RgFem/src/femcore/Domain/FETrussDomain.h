#pragma once
#include "FEBeamDomain.h"

class FETrussElement;

//-----------------------------------------------------------------------------
//! Abstract base class for truss elements
class FEM_EXPORT FETrussDomain : public FEBeamDomain
{
public:
	FETrussDomain(FEModel* pm);

public:
	bool Create(int nsize, FE_Element_Spec espec) override;

	int Elements() const override { return (int)m_Elem.size(); }

	FETrussElement& Element(int i) { return m_Elem[i]; }

	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

public:
	void ForEachTrussElement(std::function<void(FETrussElement& el)> f);

public:
	//! Calculate the truss normal
	Vector3d TrussNormal(FETrussElement& el);

protected:
	std::vector<FETrussElement>	m_Elem;
};
