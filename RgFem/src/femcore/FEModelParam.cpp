#include "FEModelParam.h"
//#include "MObjBuilder.h"
#include "FEDataArray.h"
#include "basicio/DumpStream.h"
//#include "FEConstValueVec3.h"

//---------------------------------------------------------------------------------------
FEModelParam::FEModelParam()
{ 
	m_scl = 1.0;
	m_dom = 0;
}

FEModelParam::~FEModelParam()
{
}

// serialization
void FEModelParam::Serialize(DumpStream& ar)
{
	ar & m_scl;
	if (ar.IsShallow() == false)
	{
		ar & m_dom;
	}
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = RANGO_NEW<FEScalarValuator>(nullptr, "const");
}

FEParamDouble::~FEParamDouble()
{
	delete m_val;
}

FEParamDouble::FEParamDouble(const FEParamDouble& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamDouble::operator = (double v)
{
	FEConstValue* val = RANGO_NEW<FEConstValue>(nullptr,"const");
	*val->constValue() = v;
	setValuator(val);
}

void FEParamDouble::operator = (const FEParamDouble& p)
{
	if (m_val) delete m_val;
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamDouble::setValuator(FEScalarValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

// get the valuator
FEScalarValuator* FEParamDouble::valuator()
{
	return m_val;
}

// is this a const value
bool FEParamDouble::isConst() const { return m_val->isConst(); };

// get the const value (returns 0 if param is not const)
double& FEParamDouble::constValue() { assert(isConst());  return *m_val->constValue(); }
double FEParamDouble::constValue() const { assert(isConst()); return *m_val->constValue(); }

void FEParamDouble::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

bool FEParamDouble::Init()
{
	return (m_val ? m_val->Init() : true);
}

//---------------------------------------------------------------------------------------
FEParamVec3::FEParamVec3()
{
	m_val = RANGO_NEW<FEVec3dValuator>(nullptr, "vector");
}

FEParamVec3::~FEParamVec3()
{
	delete m_val;
}

FEParamVec3::FEParamVec3(const FEParamVec3& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

bool FEParamVec3::Init()
{
	return (m_val ? m_val->Init() : true);
}

// set the value
void FEParamVec3::operator = (const Vector3d& v)
{
	FEConstValueVec3* val = RANGO_NEW<FEConstValueVec3>(nullptr,"vector");
	val->value() = v;
	setValuator(val);
}

void FEParamVec3::operator = (const FEParamVec3& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamVec3::setValuator(FEVec3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

FEVec3dValuator* FEParamVec3::valuator()
{
	return m_val;
}

void FEParamVec3::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

//==========================================================================

FEParamMat3d::FEParamMat3d()
{
	m_val = RANGO_NEW<FEMat3dValuator>(nullptr, "const");
}

FEParamMat3d::~FEParamMat3d()
{
	delete m_val;
}

FEParamMat3d::FEParamMat3d(const FEParamMat3d& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamMat3d::operator = (const Matrix3d& v)
{
	FEConstValueMat3d* val = RANGO_NEW<FEConstValueMat3d>(nullptr, "const");
	val->value() = v;
	setValuator(val);
}

void FEParamMat3d::operator = (const FEParamMat3d& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

bool FEParamMat3d::Init()
{
	return (m_val ? m_val->Init() : true);
}

// set the valuator
void FEParamMat3d::setValuator(FEMat3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

// get the valuator
FEMat3dValuator* FEParamMat3d::valuator()
{
	return m_val;
}

void FEParamMat3d::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

//==========================================================================
FEParamMat3ds::FEParamMat3ds()
{
	m_val = RANGO_NEW<FEMat3dsValuator>(nullptr, "const");
}

FEParamMat3ds::~FEParamMat3ds()
{
	delete m_val;
}

FEParamMat3ds::FEParamMat3ds(const FEParamMat3ds& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamMat3ds::operator = (const Matrix3ds& v)
{
	FEConstValueMat3ds* val = RANGO_NEW<FEConstValueMat3ds>(nullptr, "const");
	val->value() = v;
	setValuator(val);
}

void FEParamMat3ds::operator = (const FEParamMat3ds& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamMat3ds::setValuator(FEMat3dsValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

FEMat3dsValuator* FEParamMat3ds::valuator()
{
	return m_val;
}

void FEParamMat3ds::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}
