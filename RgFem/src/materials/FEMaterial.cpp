
#include "FEMaterial.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEMaterialBase::FEMaterialBase(FEModel* fem) : FEModelComponent(fem)
{
}

//-----------------------------------------------------------------------------
//! returns a pointer to a new material point object
FEMaterialPointData* FEMaterialBase::CreateMaterialPointData() { return nullptr; };

//-----------------------------------------------------------------------------
//! Update specialized material points at each iteration
void FEMaterialBase::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{

}

//=============================================================================
BEGIN_PARAM_DEFINE(FEMaterial, FEMaterialBase)
	//ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);
END_PARAM_DEFINE();

//-----------------------------------------------------------------------------
FEMaterial::FEMaterial(FEModel* fem) : FEMaterialBase(fem)
{
	m_Q = nullptr;
}

//-----------------------------------------------------------------------------
FEMaterial::~FEMaterial()
{
}

//-----------------------------------------------------------------------------
// evaluate local coordinate system at material point
Matrix3d FEMaterial::GetLocalCS(const FEMaterialPoint& mp)
{
	Matrix3d Q = (m_Q ? m_Q->operator()(mp) : Matrix3d::identity());
	FEMaterial* parent = dynamic_cast<FEMaterial*>(GetParent());
	if (parent) 
	{
		Matrix3d Qp = parent->GetLocalCS(mp);
		return Qp*Q;
	}
	else
	{
		Matrix3d A = mp.m_Q.RotationMatrix();
		return A*Q;
	}
}

//-----------------------------------------------------------------------------
// set the (local) material axis valuator
void FEMaterial::SetMaterialAxis(FEMat3dValuator* val)
{
	if (m_Q) delete m_Q;
	m_Q = val;
}

//-----------------------------------------------------------------------------
//! Initial material.
bool FEMaterial::Init()
{
	// initialize base class
	return FECoreBase::Init();
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomain(FEDomain* dom)
{
	m_domList.AddDomain(dom);
}

//-----------------------------------------------------------------------------
FEDomainParameter* FEMaterial::FindDomainParameter(const std::string& paramName)
{
	for (int i = 0; i < m_param.size(); ++i)
	{
		FEDomainParameter* pi = m_param[i];
		if (pi->name() == paramName) return pi;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomainParameter(FEDomainParameter* p)
{
	assert(p);
	m_param.push_back(p);
}

//==============================================================================
BEGIN_FECORE_CLASS(FEMaterialProperty, FEMaterialBase)
END_FECORE_CLASS();

FEMaterialProperty::FEMaterialProperty(FEModel* fem) : FEMaterialBase(fem)
{

}

//-----------------------------------------------------------------------------
// Since properties don't have local coordinate system,
// we return the parent's 
Matrix3d FEMaterialProperty::GetLocalCS(const FEMaterialPoint& mp)
{
	FEMaterialBase* parent = dynamic_cast<FEMaterialBase*>(GetParent()); assert(parent);
	return parent->GetLocalCS(mp);
}
