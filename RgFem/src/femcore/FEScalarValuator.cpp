#include "femcore/FEScalarValuator.h"
#include "materials/FEMaterialPoint.h"
#include "femcore/FEModelParam.h"
#include "femcore/FEModel.h"
#include "basicio/DumpStream.h"
#include "logger/log.h"

DEFINE_META_CLASS(FEScalarValuator, FEValuator, "");

////=============================================================================
//bool FEMathExpression::Init(const std::string& expr, FEObjectBase* pc)
//{
//	Clear();
//
//	AddVariable("X");
//	AddVariable("Y");
//	AddVariable("Z");
//	AddVariable("t");
//	bool b = Create(expr, true);
//
//	// lookup all the other variables.
//	if (Variables() > 4)
//	{
//		if (pc == nullptr) return false;
//
//		for (int i = 4; i < Variables(); ++i)
//		{
//			MVariable* vari = Variable(i);
//
//			ParamString ps(vari->Name().c_str());
//			FEParam* p = pc->FindParameter(ps);
//			if (p)
//			{
//				assert((p->type() == FE_PARAM_DOUBLE_MAPPED) || (p->type() == FE_PARAM_DOUBLE) || (p->type() == FE_PARAM_INT));
//
//				MathParam mp;
//				mp.type = 0;
//				mp.pp = p;
//				m_vars.push_back(mp);
//			}
//			else
//			{
//				// see if it's a data map
//				FEModel* fem = pc->GetFEModel();
//
//				// see if it's a global variable
//				p = fem->FindParameter(ps);
//				if (p)
//				{
//					MathParam mp;
//					mp.type = 0;
//					mp.pp = p;
//					m_vars.push_back(mp);
//				}
//				else
//				{
//					FEMesh& mesh = fem->GetMesh();
//
//					FEDataMap* map = mesh.FindDataMap(vari->Name());
//					assert(map);
//					if (map == nullptr) {
//						const char* szvar = vari->Name().c_str();
//						feLogErrorEx(fem, "Don't understand variable name \"%s\" in math expression.", szvar);
//						return false;
//					}
//					if (map->DataType() != FEDataType::FE_DOUBLE) {
//						const char* szvar = vari->Name().c_str();
//						feLogErrorEx(fem, "Variable \"%s\" is not a scalar variable.", szvar);
//						return false;
//					}
//
//					MathParam mp;
//					mp.type = 1;
//					mp.map = map;
//					m_vars.push_back(mp);
//				}
//			}
//		}
//	}
//
//	assert(b);
//	return b;
//}
//
//void FEMathExpression::operator = (const FEMathExpression& me)
//{
//	MSimpleExpression::operator=(me);
//	m_vars = me.m_vars;
//}
//
//double FEMathExpression::value(FEModel* fem, const FEMaterialPoint& pt)
//{
//	std::vector<double> var(4 + m_vars.size());
//	var[0] = pt.m_r0.x;
//	var[1] = pt.m_r0.y;
//	var[2] = pt.m_r0.z;
//	var[3] = fem->GetTime().currentTime;
//	if (m_vars.empty() == false)
//	{
//		for (int i = 0; i < (int)m_vars.size(); ++i)
//		{
//			MathParam& mp = m_vars[i];
//			if (mp.type == 0)
//			{
//				FEParam* pi = mp.pp;
//				switch (pi->type())
//				{
//				case FE_PARAM_INT: var[4 + i] = (double)pi->value<int>(); break;
//				case FE_PARAM_DOUBLE: var[4 + i] = pi->value<double>(); break;
//				case FE_PARAM_DOUBLE_MAPPED: var[4 + i] = pi->value<FEParamDouble>()(pt); break;
//				}
//			}
//			else
//			{
//				FEDataMap& map = *mp.map;
//				var[4 + i] = map.value(pt);
//			}
//		}
//	}
//	return value_s(var);
//}
//
////=============================================================================
//BEGIN_FECORE_CLASS(FEConstValue, FEScalarValuator)
//	ADD_PARAMETER(m_val, "const");
//END_FECORE_CLASS();
//
//FEScalarValuator* FEConstValue::copy()
//{
//	FEConstValue* val = fecore_alloc(FEConstValue, GetFEModel());
//	val->m_val = m_val;
//	return val;
//}
//
////=============================================================================
//
//BEGIN_FECORE_CLASS(FEMathValue, FEScalarValuator)
//	ADD_PARAMETER(m_expr, "math");
//END_FECORE_CLASS();
//
//FEMathValue::FEMathValue(FEModel* fem) : FEScalarValuator(fem)
//{
//	m_parent = nullptr;
//}
//
//void FEMathValue::setMathString(const std::string& s)
//{
//	m_expr = s;
//}
//
//bool FEMathValue::Init()
//{
//	return create();
//}
//
//void FEMathValue::Serialize(DumpStream& ar)
//{
//	FEScalarValuator::Serialize(ar);
//	if (ar.IsShallow()) return;
//
//	if (ar.IsSaving())
//	{
//		ar << m_parent;
//	}
//	else if (ar.IsLoading())
//	{
//		ar >> m_parent;
//		create(m_parent);
//	}
//}
//
//bool FEMathValue::create(FEObjectBase* pc)
//{
//	// see if this is already initialized
//	if (m_math.Variables() > 0) return true;
//
//	if (pc == nullptr)
//	{
//		// try to find the owner of this parameter
//		// First, we need the model parameter
//		FEModelParam* param = GetModelParam();
//		if (param == nullptr) return false;
//
//		// we'll need the model for this
//		FEModel* fem = GetFEModel();
//		if (fem == nullptr) return false;
//
//		// Now try to find the owner of this parameter
//		pc = fem->FindParameterOwner(param);
//	}
//	m_parent = pc;
//
//	// initialize the math expression
//	bool b = m_math.Init(m_expr, pc);
//	return b;
//}
//
//FEMathValue::~FEMathValue()
//{
//}

//FEScalarValuator* FEMathValue::copy()
//{
//	FEMathValue* newExpr = fecore_alloc(FEMathValue, GetFEModel());
//	newExpr->m_expr = m_expr;
//	newExpr->m_math = m_math;
//	return newExpr;
//}
//
//double FEMathValue::operator()(const FEMaterialPoint& pt)
//{
//	return m_math.value(GetFEModel(), pt);
//}

//---------------------------------------------------------------------------------------

FEMappedValue::FEMappedValue() : FEScalarValuator(), m_val(nullptr)
{
	m_scale = 1.0;
}

void FEMappedValue::setDataMap(FEDataMap* val)
{
	m_val = val;
}

FEDataMap* FEMappedValue::dataMap()
{
	return m_val;
}

void FEMappedValue::setScaleFactor(double s)
{
	m_scale = s;
}

double FEMappedValue::operator()(const FEMaterialPoint& pt)
{
	return m_scale*m_val->value(pt);
}

FEScalarValuator* FEMappedValue::copy()
{
	FEMappedValue* map = RANGO_NEW<FEMappedValue>(GetFEModel(), "");
	map->setDataMap(m_val);
	map->m_scale = m_scale;
	return map;
}

void FEMappedValue::Serialize(DumpStream& dmp)
{
	if (dmp.IsShallow()) return;
	dmp & m_scale;
	dmp & m_val;
}

//---------------------------------------------------------------------------------------

//FENodeMappedValue::FENodeMappedValue(FEModel* fem) : FEScalarValuator(fem), m_val(nullptr)
//{
//
//}
//
//void FENodeMappedValue::setDataMap(FENodeDataMap* val)
//{
//	m_val = val;
//}
//
//double FENodeMappedValue::operator()(const FEMaterialPoint& pt)
//{
//	return m_val->getValue(pt.m_index);
//}
//
//FEScalarValuator* FENodeMappedValue::copy()
//{
//	FENodeMappedValue* map = RANGO_NEW<FENodeMappedValue>(GetFEModel(),"");
//	map->setDataMap(m_val);
//	return map;
//}
