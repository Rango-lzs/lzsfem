#pragma once
#include "FEModelLoad.h"
#include "femcore/Domain/RgDomain.h"
#include "femcore/Domain/RgDomainList.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Base class for body-loads
class FEM_EXPORT FEBodyLoad : public FEModelLoad
{
	DECLARE_META_CLASS(FEBodyLoad, FEModelLoad);

public:
	FEBodyLoad();
	virtual ~FEBodyLoad();

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

public:
	//! return number of domains this load is applied to
	int Domains() const;

	//! return a domain 
	RgDomain* Domain(int i);

	//! add a domain to which to apply this load
	void SetDomainList(RgElementSet* elset);

	//! get the domain list
	RgDomainList& GetDomainList();
    
private:
	RgDomainList	m_dom;	//!< list of domains to which to apply the body load
};
