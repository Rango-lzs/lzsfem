#pragma once
#include "femcore/FEValuator.h"
//#include "datastructure/MathObject.h"
#include "femcore/FEDataMap.h"
#include "femcore/FENodeDataMap.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating scalar parameters
class FEM_EXPORT FEScalarValuator : public FEValuator
{
    DECLARE_META_CLASS(FEScalarValuator, FEValuator);

public:
	FEScalarValuator(FEModel* fem) : FEValuator(fem) {};

	virtual double operator()(const FEMaterialPoint& pt) = 0;

	virtual FEScalarValuator* copy() = 0;

	virtual bool isConst() { return false; }

	virtual double* constValue() { return nullptr; }
};

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEConstValue : public FEScalarValuator
{
public:
	FEConstValue(FEModel* fem) : FEScalarValuator(fem), m_val(0.0) {};
	double operator()(const FEMaterialPoint& pt) override { return m_val; }

	bool isConst() override { return true; }

	double* constValue() override { return &m_val; }

	FEScalarValuator* copy() override;

private:
	double	m_val;

	DECLARE_PARAM_LIST();
};

////---------------------------------------------------------------------------------------
//class FEMathExpression : public MSimpleExpression
//{
//	struct MathParam
//	{
//		int			type;	// 0 = param, 1 = map
//		FEParam* pp;
//		FEDataMap* map;
//	};
//
//public:
//	bool Init(const std::string& expr, FEObjectBase* pc = nullptr);
//
//	void operator = (const FEMathExpression& me);
//
//	double value(FEModel* fem, const FEMaterialPoint& pt);
//
//private:
//	std::vector<MathParam>	m_vars;
//};
//
////---------------------------------------------------------------------------------------
//class FEM_EXPORT FEMathValue : public FEScalarValuator
//{
//public:
//	FEMathValue(FEModel* fem);
//	~FEMathValue();
//	double operator()(const FEMaterialPoint& pt) override;
//
//	bool Init() override;
//
//	FEScalarValuator* copy() override;
//
//	void setMathString(const std::string& s);
//
//	bool create(FEObjectBase* pc = 0);
//
//	void Serialize(DumpStream& ar) override;
//
//private:
//	std::string			m_expr;
//	FEMathExpression	m_math;
//	FEObjectBase*			m_parent;
//
//	DECLARE_FECORE_CLASS();
//};

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEMappedValue : public FEScalarValuator
{
public:
	FEMappedValue(FEModel* fem);
	void setDataMap(FEDataMap* val);
	void setScaleFactor(double s);

	FEDataMap* dataMap();

	double operator()(const FEMaterialPoint& pt) override;

	FEScalarValuator* copy() override;

	void Serialize(DumpStream& dmp) override;

private:
	double		m_scale;	// scale factor
	FEDataMap*	m_val;
};

////---------------------------------------------------------------------------------------
//class FEM_EXPORT FENodeMappedValue : public FEScalarValuator
//{
//public:
//	FENodeMappedValue(FEModel* fem);
//	void setDataMap(FENodeDataMap* val);
//
//	double operator()(const FEMaterialPoint& pt) override;
//
//	FEScalarValuator* copy() override;
//
//private:
//	FENodeDataMap*		m_val;
//};

