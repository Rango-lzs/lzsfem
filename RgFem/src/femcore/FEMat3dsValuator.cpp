#include "FEMat3dsValuator.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "FEDataMap.h"
#include "basicio/DumpStream.h"

//=============================================================================
// FEConstValueMat3ds
//-----------------------------------------------------------------------------

BEGIN_PARAM_DEFINE(FEConstValueMat3ds, FEMat3dsValuator)
	ADD_PARAMETER(m_val, "const");
END_PARAM_DEFINE();

FEConstValueMat3ds::FEConstValueMat3ds(FEModel* fem) : FEMat3dsValuator(fem)
{
	m_val.zero();
}

FEMat3dsValuator* FEConstValueMat3ds::copy()
{
    FEConstValueMat3ds* map = RANGO_NEW<FEConstValueMat3ds>(GetFEModel(), "");
	map->m_val = m_val;
	return map;
}

//=============================================================================
// FEMappedValueMat3ds
//-----------------------------------------------------------------------------

BEGIN_PARAM_DEFINE(FEMappedValueMat3ds, FEMat3dsValuator)
	ADD_PARAMETER(m_mapName, "map");
END_PARAM_DEFINE();

FEMappedValueMat3ds::FEMappedValueMat3ds(FEModel* fem) : FEMat3dsValuator(fem)
{
	m_val = nullptr;
}

void FEMappedValueMat3ds::setDataMap(FEDataMap* val)
{
	m_val = val;
}

Matrix3ds FEMappedValueMat3ds::operator()(const FEMaterialPoint& pt)
{
	return m_val->valueMat3ds(pt);
}

FEMat3dsValuator* FEMappedValueMat3ds::copy()
{
    FEMappedValueMat3ds* map = RANGO_NEW<FEMappedValueMat3ds>(GetFEModel(), "");
	map->m_val = m_val;
	return map;
}

void FEMappedValueMat3ds::Serialize(DumpStream& ar)
{
	FEMat3dsValuator::Serialize(ar);
	ar & m_val;
}
