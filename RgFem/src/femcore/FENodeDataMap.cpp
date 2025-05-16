#include "FENodeDataMap.h"
#include "FENodeSet.h"
#include "materials/FEMaterialPoint.h"
#include "basicio/DumpStream.h"

FENodeDataMap::FENodeDataMap() : FEDataMap(FE_NODE_DATA_MAP, FE_INVALID_TYPE)
{
	m_nodeSet = nullptr;
}

FENodeDataMap::FENodeDataMap(FEDataType dataType) : FEDataMap(FE_NODE_DATA_MAP, dataType)
{
	m_nodeSet = nullptr;
}

void FENodeDataMap::Create(const FENodeSet* nodeSet, double val)
{
	m_nodeSet = nodeSet;
	int nsize = nodeSet->Size();
	resize(nsize, val);
}

const FENodeSet* FENodeDataMap::GetNodeSet() const
{ 
	return m_nodeSet; 
}

double FENodeDataMap::getValue(int n) const
{
	return get<double>(n);
}

void FENodeDataMap::setValue(int n, double v)
{
	set<double>(n, v);
}

void FENodeDataMap::setValue(int n, const Vector2d& v)
{
	set<Vector2d>(n, v);
}

void FENodeDataMap::setValue(int n, const Vector3d& v)
{
	set<Vector3d>(n, v);
}

void FENodeDataMap::setValue(int n, const Matrix3d& v)
{
	set<Matrix3d>(n, v);
}

void FENodeDataMap::setValue(int n, const Matrix3ds& v)
{
	set<Matrix3ds>(n, v);
}

void FENodeDataMap::fillValue(double v)
{
	set<double>(v);
}

void FENodeDataMap::fillValue(const Vector2d& v)
{
	set<Vector2d>(v);
}

void FENodeDataMap::fillValue(const Vector3d& v)
{
	set<Vector3d>(v);
}

void FENodeDataMap::fillValue(const Matrix3d& v)
{
	set<Matrix3d>(v);
}

void FENodeDataMap::fillValue(const Matrix3ds& v)
{
	set<Matrix3ds>(v);
}

double FENodeDataMap::value(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<double>(mp.m_index);
}

Vector3d FENodeDataMap::valueVec3d(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<Vector3d>(mp.m_index);
}

Matrix3d FENodeDataMap::valueMat3d(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<Matrix3d>(mp.m_index);
}

Matrix3ds FENodeDataMap::valueMat3ds(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<Matrix3ds>(mp.m_index);
}

// return the item list associated with this map
FEItemList* FENodeDataMap::GetItemList()
{
	return const_cast<FENodeSet*>(m_nodeSet);
}

void FENodeDataMap::Serialize(DumpStream& ar)
{
	FEDataMap::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			// We have to cast the const away before serializing
			FENodeSet* ns = const_cast<FENodeSet*>(m_nodeSet);
			ar << ns;
		}
		else
		{
			FENodeSet* ns;
			ar >> ns;
			m_nodeSet = ns;
		}
	}
}
