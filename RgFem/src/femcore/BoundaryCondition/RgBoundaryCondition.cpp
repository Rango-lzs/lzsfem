#include "RgBoundaryCondition.h"
#include "femcore/FELoadController.h"
#include "femcore/FEMesh.h"
#include "femcore/FEModel.h"
#include "femcore/FENodeSet.h"
#include "logger/log.h"

#include <algorithm>

//=============================================================================
// RgBoundaryCondition
//=============================================================================

DEFINE_META_CLASS(RgBoundaryCondition, FEObjectBase, "");

RgBoundaryCondition::RgBoundaryCondition(FEModel* fem)
    : m_fem(fem)
    , m_nodeSet(nullptr)
{

}

RgBoundaryCondition::~RgBoundaryCondition()
{
}

bool RgBoundaryCondition::Init()
{
    if (m_nodeSet == nullptr)
    {
        feLogError("Node set not set for boundary condition %s", GetName().c_str());
        return false;
    }

    return FEObjectBase::Init();
}


void RgBoundaryCondition::UpdateModel()
{
    // Default: do nothing
}

void RgBoundaryCondition::SetNodeSet(FENodeSet* nodeSet)
{
    m_nodeSet = nodeSet;
}

void RgBoundaryCondition::Serialize(DumpStream& ar)
{
    FEObjectBase::Serialize(ar);

    if (ar.IsShallow())
        return;

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
// RgFixedBC
//=============================================================================

DEFINE_META_CLASS(RgFixedBC, RgBoundaryCondition, "");

RgFixedBC::RgFixedBC(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_dofMask(DOF_ALL)
{
}

bool RgFixedBC::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    return true;
}

void RgFixedBC::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);

        // Fix DOFs based on mask
        if (m_dofMask & DOF_X)
            node.setDofState(0, DOF_FIXED);
        if (m_dofMask & DOF_Y)
            node.setDofState(1, DOF_FIXED);
        if (m_dofMask & DOF_Z)
            node.setDofState(2, DOF_FIXED);
        if (m_dofMask & DOF_RX)
            node.setDofState(3, DOF_FIXED);
        if (m_dofMask & DOF_RY)
            node.setDofState(4, DOF_FIXED);
        if (m_dofMask & DOF_RZ)
            node.setDofState(5, DOF_FIXED);

        // Set values to zero
        if (m_dofMask & DOF_X)
            node.set(0, 0.0);
        if (m_dofMask & DOF_Y)
            node.set(1, 0.0);
        if (m_dofMask & DOF_Z)
            node.set(2, 0.0);
        if (m_dofMask & DOF_RX)
            node.set(3, 0.0);
        if (m_dofMask & DOF_RY)
            node.set(4, 0.0);
        if (m_dofMask & DOF_RZ)
            node.set(5, 0.0);
    }
}

void RgFixedBC::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_dofMask;
}

//=============================================================================
// RgPrescribedDisplacement
//=============================================================================

DEFINE_META_CLASS(RgPrescribedDisplacement, RgBoundaryCondition, "");

RgPrescribedDisplacement::RgPrescribedDisplacement(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_dof(0)
    , m_scale(0.0)
    , m_loadController(nullptr)
{
}

bool RgPrescribedDisplacement::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    if (m_dof < 0 || m_dof > 2)
    {
        feLogError("Invalid DOF index %d for prescribed displacement", m_dof);
        return false;
    }

    return true;
}

void RgPrescribedDisplacement::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);

        // Mark DOF as prescribed
        node.setDofState(m_dof, DOF_PRESCRIBED);
    }
}

void RgPrescribedDisplacement::PrepStep()
{
    // Store initial values if relative
    if (!m_nodeSet)
        return;

    // Could store reference values here for relative BCs
}

void RgPrescribedDisplacement::UpdateModel()
{
    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    const FETimeInfo& tp = m_fem->GetTime();

    // Get current value from load controller
    double value = m_scale;
    if (m_loadController)
    {
        value = m_scale * m_loadController->Value();
    }

    int nnodes = m_nodeSet->Size();
    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.set(m_dof, value);
    }
}

void RgPrescribedDisplacement::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_dof;
    ar& m_scale;

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
            m_loadController = m_fem->GetLoadController(lcID);
        }*/
    }
}

//=============================================================================
// RgPrescribedRotation
//=============================================================================

DEFINE_META_CLASS(RgPrescribedRotation, RgBoundaryCondition, "");

RgPrescribedRotation::RgPrescribedRotation(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_dof(3)
    , m_scale(0.0)
    , m_loadController(nullptr)
{
}

bool RgPrescribedRotation::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    if (m_dof < 3 || m_dof > 5)
    {
        feLogError("Invalid DOF index %d for prescribed rotation", m_dof);
        return false;
    }

    return true;
}

void RgPrescribedRotation::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.setDofState(m_dof, DOF_PRESCRIBED);
    }
}

void RgPrescribedRotation::UpdateModel()
{
    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();

    double value = m_scale;
    if (m_loadController)
    {
        value = m_scale * m_loadController->Value();
    }

    int nnodes = m_nodeSet->Size();
    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.set(m_dof, value);
    }
}

void RgPrescribedRotation::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_dof;
    ar& m_scale;
}

//=============================================================================
// FEPrescribedVelocity
//=============================================================================

DEFINE_META_CLASS(RgPrescribedVelocity, RgBoundaryCondition, "");

RgPrescribedVelocity::RgPrescribedVelocity(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_dof(0)
    , m_scale(0.0)
    , m_loadController(nullptr)
{
}

bool RgPrescribedVelocity::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;
    return true;
}

void RgPrescribedVelocity::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.setDofState(m_dof, DOF_PRESCRIBED);
    }
}

void RgPrescribedVelocity::UpdateModel()
{
    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();

    double value = m_scale;
    if (m_loadController)
    {
        value = m_scale * m_loadController->Value();
    }

    int nnodes = m_nodeSet->Size();
    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        // Set velocity
        //node.m_vt = node.m_vp;  // Update velocity
    }
}

void RgPrescribedVelocity::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_dof;
    ar& m_scale;
}

//=============================================================================
// PrescribedDOF
//=============================================================================

DEFINE_META_CLASS(RgPrescribedDOF, RgBoundaryCondition, "");

RgPrescribedDOF::RgPrescribedDOF(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_dof(-1)
    , m_value(0.0)
    , m_loadController(nullptr)
{
}

bool RgPrescribedDOF::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    // Resolve DOF name to index if needed
    if (m_dof < 0 && !m_dofName.empty())
    {
        m_dof = m_fem->GetDOFIndex(m_dofName.c_str());
        if (m_dof < 0)
        {
            feLogError("Unknown DOF name: %s", m_dofName.c_str());
            return false;
        }
    }

    return true;
}

void RgPrescribedDOF::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.setDofState(m_dof, DOF_PRESCRIBED);
    }
}

void RgPrescribedDOF::PrepStep()
{
    // Could implement reference value storage for relative BCs
}

void RgPrescribedDOF::UpdateModel()
{
    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();

    double value = m_value;
    if (m_loadController)
    {
        value = m_value * m_loadController->Value();
    }

    int nnodes = m_nodeSet->Size();
    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);
        node.set(m_dof, value);
    }
}

void RgPrescribedDOF::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_dof;
    ar& m_dofName;
    ar& m_value;
}

//=============================================================================
// FESymmetryBC
//=============================================================================

DEFINE_META_CLASS(RgSymmetryBC, RgBoundaryCondition, "");

RgSymmetryBC::RgSymmetryBC()
    : RgBoundaryCondition(nullptr)
    , m_normal(0, 0, 1)
    , m_dofMask(0)
{
}

bool RgSymmetryBC::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    // Normalize the normal vector
    m_normal.unit();

    // Determine which DOFs to constrain based on normal
    m_dofMask = 0;

    // For a symmetry plane with normal n, constrain:
    // - Displacement in normal direction
    // - Rotations perpendicular to normal

    double tol = 1e-6;
    if (fabs(m_normal.x) > tol)
        m_dofMask |= RgFixedBC::DOF_X;
    if (fabs(m_normal.y) > tol)
        m_dofMask |= RgFixedBC::DOF_Y;
    if (fabs(m_normal.z) > tol)
        m_dofMask |= RgFixedBC::DOF_Z;

    return true;
}

void RgSymmetryBC::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);

        // Fix normal displacement
        if (m_dofMask & RgFixedBC::DOF_X)
        {
            node.setDofState(0, DOF_FIXED);
            node.set(0, 0.0);
        }
        if (m_dofMask & RgFixedBC::DOF_Y)
        {
            node.setDofState(1, DOF_FIXED);
            node.set(1, 0.0);
        }
        if (m_dofMask & RgFixedBC::DOF_Z)
        {
            node.setDofState(2, DOF_FIXED);
            node.set(2, 0.0);
        }
    }
}

void RgSymmetryBC::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_normal;
    ar& m_dofMask;
}

//=============================================================================
// FEAntiSymmetryBC
//=============================================================================

DEFINE_META_CLASS(RgAntiSymmetryBC, RgBoundaryCondition, "");

RgAntiSymmetryBC::RgAntiSymmetryBC(FEModel* fem)
    : RgBoundaryCondition(fem)
    , m_normal(0, 0, 1)
    , m_dofMask(0)
{
}

 RgAntiSymmetryBC::RgAntiSymmetryBC()
     :RgBoundaryCondition(nullptr)
{
}

bool RgAntiSymmetryBC::Init()
{
    if (!RgBoundaryCondition::Init())
        return false;

    m_normal.unit();

    // For anti-symmetry: constrain tangential displacements
    m_dofMask = RgFixedBC::DOF_ALL;

    double tol = 1e-6;
    if (fabs(m_normal.x) > tol)
        m_dofMask &= ~RgFixedBC::DOF_X;
    if (fabs(m_normal.y) > tol)
        m_dofMask &= ~RgFixedBC::DOF_Y;
    if (fabs(m_normal.z) > tol)
        m_dofMask &= ~RgFixedBC::DOF_Z;

    return true;
}

void RgAntiSymmetryBC::Activate()
{
    RgBoundaryCondition::Activate();

    if (!m_nodeSet)
        return;

    FEMesh& mesh = m_fem->GetMesh();
    int nnodes = m_nodeSet->Size();

    for (int i = 0; i < nnodes; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);

        if (m_dofMask & RgFixedBC::DOF_X)
        {
            node.setDofState(0, DOF_FIXED);
            node.set(0, 0.0);
        }
        if (m_dofMask & RgFixedBC::DOF_Y)
        {
            node.setDofState(1, DOF_FIXED);
            node.set(1, 0.0);
        }
        if (m_dofMask & RgFixedBC::DOF_Z)
        {
            node.setDofState(2, DOF_FIXED);
            node.set(2, 0.0);
        }
    }
}

void RgAntiSymmetryBC::Serialize(DumpStream& ar)
{
    RgBoundaryCondition::Serialize(ar);
    ar& m_normal;
    ar& m_dofMask;
}

//=============================================================================
// FEBCFactory
//=============================================================================

RgBoundaryCondition* RgBCFactory::Create(FEModel* fem, const std::string& type)
{
    if (type == "ENCASTRE" || type == "encastre" || type == "fixed")
        return new RgFixedBC(fem);
    else if (type == "displacement")
        return new RgPrescribedDisplacement(fem);
    else if (type == "rotation")
        return new RgPrescribedRotation(fem);
    else if (type == "velocity")
        return new RgPrescribedVelocity(fem);
    else if (type == "symmetry")
        return new RgSymmetryBC();
    else if (type == "antisymmetry")
        return new RgAntiSymmetryBC();
    else if (type == "prescribed")
        return new RgPrescribedDOF(fem);

    return nullptr;
}

RgFixedBC* RgBCFactory::CreateFixed(FEModel* fem, FENodeSet* nodeSet, int dofMask)
{
    RgFixedBC* bc = new RgFixedBC(fem);
    bc->SetNodeSet(nodeSet);
    bc->SetDOFs(dofMask);
    return bc;
}

RgPrescribedDisplacement* RgBCFactory::CreatePrescribed(FEModel* fem, FENodeSet* nodeSet, int dof, double value)
{
    RgPrescribedDisplacement* bc = new RgPrescribedDisplacement(fem);
    bc->SetNodeSet(nodeSet);
    bc->SetDOF(dof);
    bc->SetScale(value);
    return bc;
}

RgFixedBC* RgBCFactory::CreateEncastre(FEModel* fem, FENodeSet* nodeSet)
{
    return CreateFixed(fem, nodeSet, RgFixedBC::DOF_ALL);
}

RgSymmetryBC* RgBCFactory::CreateSymmetry(FEModel* fem, FENodeSet* nodeSet, const Vector3d& normal)
{
    RgSymmetryBC* bc = new RgSymmetryBC();
    bc->SetNodeSet(nodeSet);
    bc->SetNormal(normal);
    return bc;
}