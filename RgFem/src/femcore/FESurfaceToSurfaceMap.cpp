#include "FESurfaceToSurfaceMap.h"
#include "FEMesh.h"
#include "FESurface.h"
#include "femcore/FEDomainMap.h"

BEGIN_PARAM_DEFINE(FESurfaceToSurfaceMap, FEElemDataGenerator)
	ADD_PROPERTY(m_func, "function");
	ADD_PROPERTY(m_surf1, "bottom_surface", FEProperty::Reference);
	ADD_PROPERTY(m_surf2, "top_surface", FEProperty::Reference);
END_PARAM_DEFINE();

FESurfaceToSurfaceMap::FESurfaceToSurfaceMap(FEModel* fem) : FEElemDataGenerator(fem)
{
	m_ccp1 = 0;
	m_ccp2 = 0;
	m_func = 0;

	m_surf1 = nullptr;
	m_surf2 = nullptr;

	m_binverted = false;
}

FESurfaceToSurfaceMap::~FESurfaceToSurfaceMap()
{
	if (m_ccp1) delete m_ccp1;
	if (m_ccp2) delete m_ccp2;
}

bool FESurfaceToSurfaceMap::Init()
{
	FEMesh& mesh = GetMesh();
	if (m_func == 0) return false;
	if ((m_surf1 == nullptr) || (m_surf2 == nullptr)) return false;
	
	// we need to invert the second surface, otherwise the normal projections won't work
	if (m_binverted == false)
	{
		m_surf2->Invert();
		m_binverted = true;
	}

	// initialize projections
	if (m_ccp1 == nullptr)
	{
		m_ccp1 = new FEClosestPointProjection(*m_surf1);
		m_ccp1->HandleSpecialCases(true);
		m_ccp1->AllowBoundaryProjections(true);
		if (m_ccp1->Init() == false) return false;
	}

	if (m_ccp2 == nullptr)
	{
		m_ccp2 = new FEClosestPointProjection(*m_surf2);
		m_ccp2->HandleSpecialCases(true);
		m_ccp2->AllowBoundaryProjections(true);
		if (m_ccp2->Init() == false) return false;
	}
    
    m_func->Init();

	return FEMeshDataGenerator::Init();
}

void FESurfaceToSurfaceMap::value(const Vector3d& x, double& data)
{
	Vector3d r(x);

	// project x onto surface 1
	Vector3d q1(0,0,0), q2(0,0,0);
	Vector2d r1, r2;
	FESurfaceElement* pe1 = m_ccp1->Project(r, q1, r1);
	if (pe1 == nullptr)
	{
		assert(false);
		data = 0;
		return;
	}

	// project x onto surface 2
	FESurfaceElement* pe2 = m_ccp2->Project(r, q2, r2);
	if (pe2 == nullptr)
	{
		assert(false);
		data = 0;
		return;
	}

	double L1 = (x - q1).norm();
	double L2 = (q2 - x).norm();

	double D = L1 + L2;
	if (D == 0.0) D = 1.0;

	// find the fractional distance
	double w = L1 / D;

	// evaluate the function
	data = m_func->value(w);
}

FEDomainMap* FESurfaceToSurfaceMap::Generate()
{
	FEElementSet* elset = GetElementSet();
	if (elset == nullptr) return nullptr;

	FEDomainMap* map = new FEDomainMap(FEDataType::FE_DOUBLE, Storage_Fmt::FMT_MATPOINTS);
	map->Create(elset);
	if (FEElemDataGenerator::Generate(*map) == false)
	{
		delete map; map = nullptr;
	}
	return map;
}
