#include "RgLoad.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENodeSet.h"
#include "femcore/FEFacetSet.h"
#include "femcore/FESurface.h"
#include "femcore/RgLoadController.h"
#include "logger/log.h"
#include "materials/RgMaterialPoint.h"

//=============================================================================
// RgLoad
//=============================================================================

DEFINE_META_CLASS(RgLoad, FEStepComponent, "");

RgLoad::RgLoad(FEModel* fem)
    : m_fem(fem), m_magnitude(0.0), m_loadController(nullptr)
{
    //SetActive(true);
}

RgLoad::~RgLoad()
{
}

bool RgLoad::Init()
{
    return FEStepComponent::Init();
}

void RgLoad::Activate()
{
    FEStepComponent::Activate();
}

void RgLoad::Deactivate()
{
    FEStepComponent::Deactivate();
}

void RgLoad::Update()
{
    // Default: do nothing
}

void RgLoad::Serialize(DumpStream& ar)
{
    FEStepComponent::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_magnitude;
    
    // Serialize load controller reference
    if (ar.IsSaving())
    {
        int lcID = m_loadController ? m_loadController->GetID() : -1;
        ar << lcID;
    }
    else
    {
        int lcID;
        ar >> lcID;
        /*if (lcID >= 0 && lcID < m_fem->LoadControllers())
        {
            m_loadController = dynamic_cast<RgLoadController*>(m_fem->GetLoadController(lcID));
        }*/
    }
}

//=============================================================================
// RgNodalLoad
//=============================================================================

DEFINE_META_CLASS(RgNodalLoad, RgLoad, "");

RgNodalLoad::RgNodalLoad(FEModel* fem)
    : RgLoad(fem), m_nodeSet(nullptr), m_force(0, 0, 0), m_dof(-1)
{
}

bool RgNodalLoad::Init()
{
    if (!RgLoad::Init()) return false;
    
    if (!m_nodeSet)
    {
        feLogError("Node set not set for nodal load %s", GetName().c_str());
        return false;
    }
    
    return true;
}

void RgNodalLoad::Activate()
{
    RgLoad::Activate();
}

void RgNodalLoad::Update()
{
    if (!m_nodeSet) return;
    
    FEMesh& mesh = m_fem->GetMesh();
    
    // Get current magnitude from load controller
    double scale = m_magnitude;
    if (m_loadController)
    {
        scale = m_magnitude * m_loadController->Value();
    }
    
    int nnodes = m_nodeSet->Size();
    
    if (m_dof >= 0)
    {
        // Single DOF load
        for (int i = 0; i < nnodes; ++i)
        {
            FENode& node = *m_nodeSet->Node(i);
            node.set_load(m_dof, scale);
        }
    }
    else
    {
        // Vector load
        Vector3d F = m_force * scale;
        for (int i = 0; i < nnodes; ++i)
        {
            FENode& node = *m_nodeSet->Node(i);
            node.set_load(0, F.x);
            node.set_load(1, F.y);
            node.set_load(2, F.z);
        }
    }
}

void RgNodalLoad::Serialize(DumpStream& ar)
{
    RgLoad::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_force;
    ar & m_dof;
    
    // Serialize node set reference
    if (ar.IsSaving())
    {
        std::string setName = m_nodeSet ? m_nodeSet->GetName() : "";
        ar << setName;
    }
    else
    {
        std::string setName;
        ar >> setName;
        if (!setName.empty())
        {
            m_nodeSet = m_fem->GetMesh().FindNodeSet(setName);
        }
    }
}

//=============================================================================
// RgSurfaceLoad
//=============================================================================

DEFINE_META_CLASS(RgSurfaceLoad, RgLoad, "");

RgSurfaceLoad::RgSurfaceLoad(FEModel* fem)
    : RgLoad(fem), m_surface(nullptr), m_facetSet(nullptr),
      m_loadType(PRESSURE), m_traction(0, 0, 0), m_bFollower(false)
{
}

bool RgSurfaceLoad::Init()
{
    if (!RgLoad::Init()) return false;
    
    if (!m_surface && !m_facetSet)
    {
        feLogError("Surface/Facet set not set for surface load %s", GetName().c_str());
        return false;
    }
    
    // Create surface from facet set if needed
    if (!m_surface && m_facetSet)
    {
        m_surface = m_fem->GetMesh().CreateSurface(*m_facetSet);
    }
    
    return true;
}

void RgSurfaceLoad::Activate()
{
    RgLoad::Activate();
}

void RgSurfaceLoad::Update()
{
    // Surface loads are typically applied during element assembly
    // This just updates the current magnitude
}

void RgSurfaceLoad::Serialize(DumpStream& ar)
{
    RgLoad::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    //ar & m_loadType;
    ar & m_traction;
    ar & m_bFollower;
    
    // Serialize surface/facet set references
    if (ar.IsSaving())
    {
        std::string facetName = m_facetSet ? m_facetSet->GetName() : "";
        ar << facetName;
    }
    else
    {
        std::string facetName;
        ar >> facetName;
        if (!facetName.empty())
        {
            m_facetSet = m_fem->GetMesh().FindFacetSet(facetName);
        }
    }
}

//=============================================================================
// RgBodyLoad
//=============================================================================

DEFINE_META_CLASS(RgBodyLoad, RgLoad, "");

RgBodyLoad::RgBodyLoad(FEModel* fem)
    : RgLoad(fem), m_bodyLoadType(GRAVITY), m_force(0, 0, 0),
      m_axis(0, 0, 1), m_origin(0, 0, 0), m_omega(0.0)
{
}

bool RgBodyLoad::Init()
{
    if (!RgLoad::Init()) return false;
    
    // Normalize axis for centrifugal load
    if (m_bodyLoadType == CENTRIFUGAL)
    {
        m_axis.unit();
    }
    
    return true;
}

void RgBodyLoad::Activate()
{
    RgLoad::Activate();
}

void RgBodyLoad::Update()
{
    // Body loads are typically applied during element assembly
}

Vector3d RgBodyLoad::EvaluateForce(const RgMaterialPoint& mp) const
{
    // Get scale from load controller
    double scale = m_magnitude;
    if (m_loadController)
    {
        scale = m_magnitude * m_loadController->Value();
    }
    
    switch (m_bodyLoadType)
    {
        case GRAVITY:
        case CONSTANT:
            return m_force * scale;
            
        case CENTRIFUGAL:
        {
            // F = rho * omega^2 * r
            // where r is position vector from axis
            Vector3d pos = mp.m_rt;  // Current position
            
            // Project position onto axis to get point on axis
            double d = (pos - m_origin) * m_axis;
            Vector3d pointOnAxis = m_origin + m_axis * d;
            
            // Radial vector from axis
            Vector3d r = pos - pointOnAxis;
            double radius = r.norm();
            
            if (radius < 1e-12)
                return Vector3d(0, 0, 0);
            
            // Centrifugal force (outward)
            double omega = m_omega * scale;
            Vector3d F = r * (omega * omega / radius);
            
            return F;
        }
        
        default:
            return Vector3d(0, 0, 0);
    }
}

void RgBodyLoad::Serialize(DumpStream& ar)
{
    RgLoad::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    //ar & m_bodyLoadType;
    ar & m_force;
    ar & m_axis;
    ar & m_origin;
    ar & m_omega;
}

//=============================================================================
// RgMomentLoad
//=============================================================================

DEFINE_META_CLASS(RgMomentLoad, RgLoad, "");

RgMomentLoad::RgMomentLoad(FEModel* fem)
    : RgLoad(fem), m_nodeSet(nullptr), m_moment(0, 0, 0)
{
}

bool RgMomentLoad::Init()
{
    if (!RgLoad::Init()) return false;
    
    if (!m_nodeSet)
    {
        feLogError("Node set not set for moment load %s", GetName().c_str());
        return false;
    }
    
    return true;
}

void RgMomentLoad::Activate()
{
    RgLoad::Activate();
}

void RgMomentLoad::Update()
{
    if (!m_nodeSet) return;
    
    FEMesh& mesh = m_fem->GetMesh();
    
    // Get current magnitude from load controller
    double scale = m_magnitude;
    if (m_loadController)
    {
        scale = m_magnitude * m_loadController->Value();
    }
    
    Vector3d M = m_moment * scale;
    
    int nnodes = m_nodeSet->Size();
    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        
        // Apply moment (DOF 3, 4, 5)
        node.set_load(3, M.x);
        node.set_load(4, M.y);
        node.set_load(5, M.z);
    }
}

void RgMomentLoad::Serialize(DumpStream& ar)
{
    RgLoad::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_moment;
    
    // Serialize node set reference
    if (ar.IsSaving())
    {
        std::string setName = m_nodeSet ? m_nodeSet->GetName() : "";
        ar << setName;
    }
    else
    {
        std::string setName;
        ar >> setName;
        if (!setName.empty())
        {
            m_nodeSet = m_fem->GetMesh().FindNodeSet(setName);
        }
    }
}

//=============================================================================
// RgLoadFactory
//=============================================================================

RgLoad* RgLoadFactory::Create(FEModel* fem, const std::string& type)
{
    if (type == "nodal" || type == "cload" || type == "concentrated")
        return new RgNodalLoad(fem);
    else if (type == "pressure" || type == "dload")
        return new RgSurfaceLoad(fem);
    else if (type == "traction")
        return new RgSurfaceLoad(fem);
    else if (type == "gravity" || type == "dload_grav")
        return new RgBodyLoad(fem);
    else if (type == "centrifugal" || type == "dload_centrif")
        return new RgBodyLoad(fem);
    else if (type == "moment")
        return new RgMomentLoad(fem);
    
    return nullptr;
}

RgNodalLoad* RgLoadFactory::CreateNodalLoad(FEModel* fem, FENodeSet* nodeSet, 
                                            const Vector3d& force)
{
    RgNodalLoad* load = new RgNodalLoad(fem);
    load->SetNodeSet(nodeSet);
    load->SetForce(force);
    load->SetMagnitude(1.0);
    return load;
}

RgNodalLoad* RgLoadFactory::CreateNodalLoad(FEModel* fem, FENodeSet* nodeSet, 
                                            int dof, double magnitude)
{
    RgNodalLoad* load = new RgNodalLoad(fem);
    load->SetNodeSet(nodeSet);
    load->SetDOF(dof);
    load->SetMagnitude(magnitude);
    return load;
}

RgSurfaceLoad* RgLoadFactory::CreatePressure(FEModel* fem, FESurface* surface, 
                                              double pressure)
{
    RgSurfaceLoad* load = new RgSurfaceLoad(fem);
    load->SetSurface(surface);
    load->SetLoadType(RgSurfaceLoad::PRESSURE);
    load->SetMagnitude(pressure);
    return load;
}

RgSurfaceLoad* RgLoadFactory::CreatePressure(FEModel* fem, FEFacetSet* facetSet, 
                                              double pressure)
{
    RgSurfaceLoad* load = new RgSurfaceLoad(fem);
    load->SetFacetSet(facetSet);
    load->SetLoadType(RgSurfaceLoad::PRESSURE);
    load->SetMagnitude(pressure);
    return load;
}

RgSurfaceLoad* RgLoadFactory::CreateTraction(FEModel* fem, FESurface* surface, 
                                              const Vector3d& traction)
{
    RgSurfaceLoad* load = new RgSurfaceLoad(fem);
    load->SetSurface(surface);
    load->SetLoadType(RgSurfaceLoad::GENERAL_TRACTION);
    load->SetTraction(traction);
    load->SetMagnitude(1.0);
    return load;
}

RgBodyLoad* RgLoadFactory::CreateGravity(FEModel* fem, const Vector3d& g)
{
    RgBodyLoad* load = new RgBodyLoad(fem);
    load->SetBodyLoadType(RgBodyLoad::GRAVITY);
    load->SetForce(g);
    load->SetMagnitude(1.0);
    return load;
}

RgBodyLoad* RgLoadFactory::CreateCentrifugal(FEModel* fem, const Vector3d& axis, 
                                              const Vector3d& origin, double omega)
{
    RgBodyLoad* load = new RgBodyLoad(fem);
    load->SetBodyLoadType(RgBodyLoad::CENTRIFUGAL);
    load->SetAxis(axis);
    load->SetOrigin(origin);
    load->SetAngularVelocity(omega);
    load->SetMagnitude(1.0);
    return load;
}

RgMomentLoad* RgLoadFactory::CreateMoment(FEModel* fem, FENodeSet* nodeSet, 
                                           const Vector3d& moment)
{
    RgMomentLoad* load = new RgMomentLoad(fem);
    load->SetNodeSet(nodeSet);
    load->SetMoment(moment);
    load->SetMagnitude(1.0);
    return load;
}
