
#include "FEMaterial.h"
#include "basicio/DumpStream.h"

DEFINE_META_CLASS(FEMaterial, FEMaterialBase, "");
DEFINE_META_CLASS(FEMaterialBase, FEModelComponent, "");

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
	return (m_Q ? m_Q->operator()(mp) : Matrix3d::identity());
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
	return FEObjectBase::Init();
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomain(FEDomain* dom)
{
	m_domList.AddDomain(dom);
}


//==============================================================================
BEGIN_PARAM_DEFINE(FEMaterialProperty, FEMaterialBase)
END_PARAM_DEFINE();

FEMaterialProperty::FEMaterialProperty(FEModel* fem) : FEMaterialBase(fem)
{

}

//-----------------------------------------------------------------------------
// Since properties don't have local coordinate system,
// we return the parent's 
Matrix3d FEMaterialProperty::GetLocalCS(const FEMaterialPoint& mp)
{
	/*FEMaterialBase* parent = dynamic_cast<FEMaterialBase*>(GetParent()); assert(parent);
	return parent->GetLocalCS(mp);*/
    return Matrix3d::identity();
}
