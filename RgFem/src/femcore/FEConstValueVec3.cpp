#include "femcore/FEModel.h"
#include "FEConstValueVec3.h"
#include "materials/FEMaterialPoint.h"
#include "FEMeshPartition.h"
#include "FENode.h"
#include "datastructure/quatd.h"
#include <assert.h>
#include "basicio/DumpStream.h"
#include "datastructure/MathUtils.h"
#include "FEMesh.h"

//==================================================================================
BEGIN_PARAM_DEFINE(FEConstValueVec3, FEVec3dValuator)
	ADD_PARAMETER(m_val, "vector");
END_PARAM_DEFINE();

FEConstValueVec3::FEConstValueVec3() : FEVec3dValuator() {}

FEVec3dValuator* FEConstValueVec3::copy()
{
	FEConstValueVec3* val = RANGO_NEW<FEConstValueVec3>( GetFEModel(),"");
	val->m_val = m_val;
	return val;
}

//==================================================================================
BEGIN_PARAM_DEFINE(FEMathValueVec3, FEVec3dValuator)
	ADD_PARAMETER(m_expr, "math");
END_PARAM_DEFINE();

FEMathValueVec3::FEMathValueVec3() : FEVec3dValuator()
{
	m_expr = "0,0,0";
	Init();
}

//---------------------------------------------------------------------------------------
bool FEMathValueVec3::Init()
{
	size_t c1 = m_expr.find(',', 0); if (c1 == std::string::npos) return false;
	size_t c2 = m_expr.find(',', c1 + 1); if (c2 == std::string::npos) return false;

	std::string sx = m_expr.substr(0, c1);
	std::string sy = m_expr.substr(c1 + 1, c2 - c1);
	std::string sz = m_expr.substr(c2 + 1, std::string::npos);

	return create(sx, sy, sz);
}

//---------------------------------------------------------------------------------------
bool FEMathValueVec3::create(const std::string& sx, const std::string& sy, const std::string& sz)
{
	FEObjectBase* pc = nullptr;
	if (pc == nullptr)
	{
		// try to find the owner of this parameter
		// First, we need the model parameter
		FEModelParam* param = GetModelParam();
		if (param == nullptr) return false;

		// we'll need the model for this
		FEModel* fem = GetFEModel();
		if (fem == nullptr) return false;

		// Now try to find the owner of this parameter
		pc = fem->FindParameterOwner(param);
	}

	/*if (m_math[0].Init(sx, pc) == false) return false;
	if (m_math[1].Init(sy, pc) == false) return false;
	if (m_math[2].Init(sz, pc) == false) return false;*/

	return true;
}

bool FEMathValueVec3::UpdateParams()
{
	return Init();
}

Vector3d FEMathValueVec3::operator()(const FEMaterialPoint& pt)
{
	/*double vx = m_math[0].value(GetFEModel(), pt);
	double vy = m_math[1].value(GetFEModel(), pt);
	double vz = m_math[2].value(GetFEModel(), pt);*/
	return Vector3d(0,0,0);
}

//---------------------------------------------------------------------------------------
FEVec3dValuator* FEMathValueVec3::copy()
{
	FEMathValueVec3* newVal = RANGO_NEW<FEMathValueVec3>( GetFEModel(),"");
	/*newVal->m_math[0] = m_math[0];
	newVal->m_math[1] = m_math[1];
	newVal->m_math[2] = m_math[2];*/
	return newVal;
}

//---------------------------------------------------------------------------------------
BEGIN_PARAM_DEFINE(FEMappedValueVec3, FEVec3dValuator)
	ADD_PARAMETER(m_mapName, "map");
END_PARAM_DEFINE();

FEMappedValueVec3::FEMappedValueVec3() : FEVec3dValuator()
{
	m_val = nullptr;
}

void FEMappedValueVec3::setDataMap(FEDataMap* val, Vector3d scl)
{
	m_val = val;
}

Vector3d FEMappedValueVec3::operator()(const FEMaterialPoint& pt)
{
	Vector3d r = m_val->valueVec3d(pt);
	return Vector3d(r.x, r.y, r.z);
}

FEVec3dValuator* FEMappedValueVec3::copy()
{
	FEMappedValueVec3* map = RANGO_NEW<FEMappedValueVec3>( GetFEModel(),"");
	map->m_val = m_val;
	return map;
}

void FEMappedValueVec3::Serialize(DumpStream& ar)
{
	FEVec3dValuator::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_val;
}

bool FEMappedValueVec3::Init()
{
	if (m_val == nullptr)
	{
		FEModel& fem = *GetFEModel();
		FEMesh& mesh = fem.GetMesh();
		FEDataMap* map = mesh.FindDataMap(m_mapName);
		if (map == nullptr) return false;
		setDataMap(map);
	}
	return FEVec3dValuator::Init();
}

//=================================================================================================
BEGIN_PARAM_DEFINE(FELocalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_n, 2, "local");
END_PARAM_DEFINE();

FELocalVectorGenerator::FELocalVectorGenerator() : FEVec3dValuator()
{
	m_n[0] = m_n[1] = 0;
}

bool FELocalVectorGenerator::Init()
{
	if ((m_n[0] <= 0) && (m_n[1] <= 0))
	{
		m_n[0] = 1;
		m_n[1] = 2;
	}
	if ((m_n[0] <= 0) || (m_n[1] <= 0)) return false;
	return FEVec3dValuator::Init();
}

Vector3d FELocalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	FEElement* el = mp.m_elem; assert(el);

	FEMeshPartition* dom = el->GetMeshPartition();
	Vector3d r0 = dom->Node(el->getLocNodeId(m_n[0]-1)).m_r0;
	Vector3d r1 = dom->Node(el->getLocNodeId(m_n[1]-1)).m_r0;

	Vector3d n = r1 - r0;
	n.unit();

	return n;
}

FEVec3dValuator* FELocalVectorGenerator::copy()
{
	FELocalVectorGenerator* map = RANGO_NEW<FELocalVectorGenerator>( GetFEModel(),"");
	map->m_n[0] = m_n[0];
	map->m_n[1] = m_n[1];
	return map;
}

//=================================================================================================
BEGIN_PARAM_DEFINE(FESphericalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_vector, "vector");
END_PARAM_DEFINE();

FESphericalVectorGenerator::FESphericalVectorGenerator() : FEVec3dValuator()
{
	m_center = Vector3d(0, 0, 0);
	m_vector = Vector3d(1, 0, 0);
}

bool FESphericalVectorGenerator::Init()
{
	// Make sure the vector is a unit vector
	m_vector.unit();
	return true;
}

FEVec3dValuator* FESphericalVectorGenerator::copy()
{
	FESphericalVectorGenerator* map = RANGO_NEW<FESphericalVectorGenerator>( GetFEModel(),"");
	map->m_center = m_center;
	map->m_vector = m_vector;
	return map;
}

Vector3d FESphericalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	Vector3d a = mp.m_r0 - m_center;
	a.unit();

	// setup the rotation
	Vector3d e1(1, 0, 0);
	quatd q(e1, a);

	Vector3d v = m_vector;
	//	v.unit();	
	q.RotateVector(v);

	return v;
}

//=================================================================================================
BEGIN_PARAM_DEFINE(FECylindricalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_axis, "axis");
	ADD_PARAMETER(m_vector, "vector");
END_PARAM_DEFINE();

FECylindricalVectorGenerator::FECylindricalVectorGenerator() : FEVec3dValuator()
{
	m_center = Vector3d(0, 0, 0);
	m_axis = Vector3d(0, 0, 1);
	m_vector = Vector3d(1, 0, 0);
}

bool FECylindricalVectorGenerator::Init()
{
	// Make sure the axis and vector are unit vectors
	m_axis.unit();
	m_vector.unit();
	return true;
}

Vector3d FECylindricalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	Vector3d p = mp.m_r0 - m_center;

	// find the vector to the axis
	Vector3d b = p - m_axis * (m_axis*p);
	b.unit();

	// setup the rotation
	Vector3d e1(1, 0, 0);
	quatd q(e1, b);

	Vector3d r = m_vector;
	//	r.unit();	
	q.RotateVector(r);

	return r;
}

FEVec3dValuator* FECylindricalVectorGenerator::copy()
{
	FECylindricalVectorGenerator* map = RANGO_NEW<FECylindricalVectorGenerator>(GetFEModel(),"");
	map->m_center = m_center;
	map->m_axis = m_axis;
	map->m_vector = m_vector;
	return map;
}


//=================================================================================================
BEGIN_PARAM_DEFINE(FESphericalAnglesVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_theta, "theta");
	ADD_PARAMETER(m_phi, "phi");
END_PARAM_DEFINE();

FESphericalAnglesVectorGenerator::FESphericalAnglesVectorGenerator() : FEVec3dValuator()
{
	// equal to x-axis (1,0,0)
	m_theta = 0.0;
	m_phi = 90.0;
}

Vector3d FESphericalAnglesVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	// convert from degress to radians
	const double the = m_theta(mp)* PI / 180.;
	const double phi = m_phi(mp)* PI / 180.;

	// the fiber vector
	Vector3d a;
	a.x = cos(the)*sin(phi);
	a.y = sin(the)*sin(phi);
	a.z = cos(phi);

	return a;
}

FEVec3dValuator* FESphericalAnglesVectorGenerator::copy()
{
	FESphericalAnglesVectorGenerator* v = RANGO_NEW<FESphericalAnglesVectorGenerator>(GetFEModel(),"");
	v->m_theta = m_theta;
	v->m_phi = m_phi;
	return v;
}


//=================================================================================================
BEGIN_PARAM_DEFINE(FEUserVectorGenerator, FEVec3dValuator)
END_PARAM_DEFINE();

FEUserVectorGenerator::FEUserVectorGenerator() : FEVec3dValuator()
{
}

Vector3d FEUserVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	assert(false);
	return Vector3d(0, 0, 0);
}

FEVec3dValuator* FEUserVectorGenerator::copy()
{
	assert(false);
	return RANGO_NEW<FEUserVectorGenerator>(GetFEModel(),"");
}
