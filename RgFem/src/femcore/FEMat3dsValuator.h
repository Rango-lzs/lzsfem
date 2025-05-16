#pragma once
#include "femcore/FEValuator.h"

class FEDataMap;
class FEMaterialPoint;

//---------------------------------------------------------------------------------------
// Base class for evaluating Vector3d parameters
class FEM_EXPORT FEMat3dsValuator : public FEValuator
{
    DECLARE_META_CLASS(FEMat3dsValuator, FEValuator);

public:
	FEMat3dsValuator(FEModel* fem) : FEValuator(fem) {};

	virtual Matrix3ds operator()(const FEMaterialPoint& pt) = 0;

	virtual FEMat3dsValuator* copy() = 0;

	virtual bool isConst() { return false; }

	virtual Matrix3ds* constValue() { return nullptr; }
};

//-----------------------------------------------------------------------------
// A constant valuator
class FEM_EXPORT FEConstValueMat3ds : public FEMat3dsValuator
{
public:
	FEConstValueMat3ds(FEModel* fem);

	FEMat3dsValuator* copy() override;

	Matrix3ds operator()(const FEMaterialPoint& pt) override { return m_val; }

	// is this a const value
	bool isConst() override { return true; }

	// get the const value (returns 0 if param is not const)
	Matrix3ds* constValue() override { return &m_val; }

	Matrix3ds& value() { return m_val; }

private:
	Matrix3ds	m_val;

	DECLARE_PARAM_LIST();
};

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEMappedValueMat3ds : public FEMat3dsValuator
{
public:
	FEMappedValueMat3ds(FEModel* fem);

	void setDataMap(FEDataMap* val);

	Matrix3ds operator()(const FEMaterialPoint& pt) override;

	FEMat3dsValuator* copy() override;

	void Serialize(DumpStream& ar) override;

private:
	std::string	m_mapName;

private:
	FEDataMap*	m_val;

	DECLARE_PARAM_LIST();
};
