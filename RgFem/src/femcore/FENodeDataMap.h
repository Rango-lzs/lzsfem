#pragma once
#include "FEDataMap.h"

class FENodeSet;

class FEM_EXPORT FENodeDataMap : public FEDataMap  
{
public:
	FENodeDataMap();
	FENodeDataMap(FEDataType dataType);

	void Create(const FENodeSet* nodeSet, double val = 0.0);

	const FENodeSet* GetNodeSet() const;

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	void Serialize(DumpStream& ar) override;

public:
	void setValue(int n, double v) override;
	void setValue(int n, const Vector2d& v) override;
	void setValue(int n, const Vector3d& v) override;
	void setValue(int n, const Matrix3d& v) override;
	void setValue(int n, const Matrix3ds& v) override;

	double getValue(int n) const;

	void fillValue(double v) override;
	void fillValue(const Vector2d& v) override;
	void fillValue(const Vector3d& v) override;
	void fillValue(const Matrix3d& v) override;
	void fillValue(const Matrix3ds& v) override;

	double value(const FEMaterialPoint& mp) override;
	Vector3d valueVec3d(const FEMaterialPoint& mp) override;
	Matrix3d valueMat3d(const FEMaterialPoint& mp) override;
	Matrix3ds valueMat3ds(const FEMaterialPoint& mp) override;

private:
	const FENodeSet*	m_nodeSet;
};
