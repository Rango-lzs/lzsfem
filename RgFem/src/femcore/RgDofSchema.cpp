/*********************************************************************
 * \file   RgDofSchema.cpp
 * \brief  Implementation of DOF Schema System
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#include "RgDofSchema.h"
#include "logger/log.h"
#include "basicio/DumpStream.h"

#include <algorithm>
#include <cassert>

 //=============================================================================
 // RgDofSet Implementation
 //=============================================================================

void RgDofSet::AddDof(DofType dof)
{
    // Check if already exists
    if (std::find(m_dofs.begin(), m_dofs.end(), dof) == m_dofs.end())
    {
        m_dofs.push_back(dof);
    }
}

bool RgDofSet::HasDof(DofType dof) const
{
    return std::find(m_dofs.begin(), m_dofs.end(), dof) != m_dofs.end();
}

DofType RgDofSet::GetDof(int index) const
{
    assert(index >= 0 && index < (int)m_dofs.size());
    return m_dofs[index];
}

//-----------------------------------------------------------------------------
// Predefined DOF Sets
//-----------------------------------------------------------------------------

RgDofSet RgDofSet::Solid3D()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    set.AddDof(DofType::DISPLACEMENT_Z);
    return set;
}

RgDofSet RgDofSet::Shell()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    set.AddDof(DofType::DISPLACEMENT_Z);
    set.AddDof(DofType::ROTATION_X);
    set.AddDof(DofType::ROTATION_Y);
    set.AddDof(DofType::ROTATION_Z);
    return set;
}

RgDofSet RgDofSet::Beam3D()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    set.AddDof(DofType::DISPLACEMENT_Z);
    set.AddDof(DofType::ROTATION_X);
    set.AddDof(DofType::ROTATION_Y);
    set.AddDof(DofType::ROTATION_Z);
    return set;
}

RgDofSet RgDofSet::Truss3D()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    set.AddDof(DofType::DISPLACEMENT_Z);
    return set;
}

RgDofSet RgDofSet::Plane2D()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    return set;
}

RgDofSet RgDofSet::Thermal3D()
{
    RgDofSet set;
    set.AddDof(DofType::TEMPERATURE);
    return set;
}

RgDofSet RgDofSet::Coupled_ThermoMech()
{
    RgDofSet set;
    set.AddDof(DofType::DISPLACEMENT_X);
    set.AddDof(DofType::DISPLACEMENT_Y);
    set.AddDof(DofType::DISPLACEMENT_Z);
    set.AddDof(DofType::TEMPERATURE);
    return set;
}

//=============================================================================
// RgDofSchema Implementation
//=============================================================================

RgDofSchema::RgDofSchema()
{
    CreateDefaultDofInfo();
}

RgDofSchema::~RgDofSchema()
{
}

//-----------------------------------------------------------------------------
void RgDofSchema::CreateDefaultDofInfo()
{
    // Create default information for all DOF types
    // This is a lookup table that can be used even if DOFs aren't active

    // Structural DOFs
    m_dofInfo[DofType::DISPLACEMENT_X] = RgDofInfo(
        DofType::DISPLACEMENT_X, DofVar::DISPLACEMENT_X, "X displacement");
    m_dofInfo[DofType::DISPLACEMENT_Y] = RgDofInfo(
        DofType::DISPLACEMENT_Y, DofVar::DISPLACEMENT_Y, "Y displacement");
    m_dofInfo[DofType::DISPLACEMENT_Z] = RgDofInfo(
        DofType::DISPLACEMENT_Z, DofVar::DISPLACEMENT_Z, "Z displacement");
    m_dofInfo[DofType::ROTATION_X] = RgDofInfo(
        DofType::ROTATION_X, DofVar::ROTATION_X, "Rotation about X");
    m_dofInfo[DofType::ROTATION_Y] = RgDofInfo(
        DofType::ROTATION_Y, DofVar::ROTATION_Y, "Rotation about Y");
    m_dofInfo[DofType::ROTATION_Z] = RgDofInfo(
        DofType::ROTATION_Z, DofVar::ROTATION_Z, "Rotation about Z");

    // Thermal DOFs
    m_dofInfo[DofType::TEMPERATURE] = RgDofInfo(
        DofType::TEMPERATURE, DofVar::TEMPERATURE, "Temperature");

    // Fluid DOFs
    m_dofInfo[DofType::VELOCITY_X] = RgDofInfo(
        DofType::VELOCITY_X, DofVar::VELOCITY_X, "X velocity");
    m_dofInfo[DofType::VELOCITY_Y] = RgDofInfo(
        DofType::VELOCITY_Y, DofVar::VELOCITY_Y, "Y velocity");
    m_dofInfo[DofType::VELOCITY_Z] = RgDofInfo(
        DofType::VELOCITY_Z, DofVar::VELOCITY_Z, "Z velocity");
    m_dofInfo[DofType::PRESSURE] = RgDofInfo(
        DofType::PRESSURE, DofVar::PRESSURE, "Pressure");

    // Electromagnetic DOFs
    m_dofInfo[DofType::ELECTRIC_POTENTIAL] = RgDofInfo(
        DofType::ELECTRIC_POTENTIAL, DofVar::ELECTRIC_POTENTIAL, "Electric potential");
    m_dofInfo[DofType::MAGNETIC_POTENTIAL_X] = RgDofInfo(
        DofType::MAGNETIC_POTENTIAL_X, DofVar::MAGNETIC_POTENTIAL_X, "Magnetic potential X");
    m_dofInfo[DofType::MAGNETIC_POTENTIAL_Y] = RgDofInfo(
        DofType::MAGNETIC_POTENTIAL_Y, DofVar::MAGNETIC_POTENTIAL_Y, "Magnetic potential Y");
    m_dofInfo[DofType::MAGNETIC_POTENTIAL_Z] = RgDofInfo(
        DofType::MAGNETIC_POTENTIAL_Z, DofVar::MAGNETIC_POTENTIAL_Z, "Magnetic potential Z");
}

//-----------------------------------------------------------------------------
void RgDofSchema::SetSchema(const RgDofSet& dofSet)
{
    Clear();

    const std::vector<DofType>& dofs = dofSet.GetDofs();
    for (DofType dof : dofs)
    {
        AddDof(dof);
    }

    RebuildIndexMaps();
}

//-----------------------------------------------------------------------------
bool RgDofSchema::AddDof(DofType dofType, const std::string& name,
    const std::string& description)
{
    // Check if already active
    if (HasDof(dofType))
    {
        RgLogWarning("DOF type %d already active", (int)dofType);
        return false;
    }

    // Get default info or create new
    RgDofInfo* info = nullptr;
    auto it = m_dofInfo.find(dofType);
    if (it != m_dofInfo.end())
    {
        info = &it->second;

        // Override with custom name/description if provided
        if (!name.empty())
            info->name = name;
        if (!description.empty())
            info->description = description;
    }
    else
    {
        // Custom DOF type - create new info
        m_dofInfo[dofType] = RgDofInfo(dofType, name, description);
        info = &m_dofInfo[dofType];
    }

    // Activate the DOF
    info->isActive = true;
    m_activeDofs.push_back(dofType);

    RebuildIndexMaps();

    return true;
}

//-----------------------------------------------------------------------------
bool RgDofSchema::AddDofs(const std::vector<DofType>& dofTypes)
{
    for (DofType dof : dofTypes)
    {
        if (!AddDof(dof))
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
bool RgDofSchema::RemoveDof(DofType dofType)
{
    if (!HasDof(dofType))
        return false;

    // Remove from active list
    auto it = std::find(m_activeDofs.begin(), m_activeDofs.end(), dofType);
    if (it != m_activeDofs.end())
    {
        m_activeDofs.erase(it);
    }

    // Mark as inactive
    auto infoIt = m_dofInfo.find(dofType);
    if (infoIt != m_dofInfo.end())
    {
        infoIt->second.isActive = false;
        infoIt->second.systemIndex = -1;
    }

    RebuildIndexMaps();

    return true;
}

//-----------------------------------------------------------------------------
void RgDofSchema::Clear()
{
    m_activeDofs.clear();
    m_nameToType.clear();
    m_typeToIndex.clear();

    // Reset all DOFs to inactive
    for (auto& pair : m_dofInfo)
    {
        pair.second.isActive = false;
        pair.second.systemIndex = -1;
    }
}

//-----------------------------------------------------------------------------
void RgDofSchema::RebuildIndexMaps()
{
    m_nameToType.clear();
    m_typeToIndex.clear();

    // Assign system indices and build maps
    for (int i = 0; i < (int)m_activeDofs.size(); ++i)
    {
        DofType dof = m_activeDofs[i];

        // Update system index in info
        auto it = m_dofInfo.find(dof);
        if (it != m_dofInfo.end())
        {
            it->second.systemIndex = i;

            // Build name mapping
            m_nameToType[it->second.name] = dof;
        }

        // Build type-to-index mapping
        m_typeToIndex[dof] = i;
    }
}

//-----------------------------------------------------------------------------
bool RgDofSchema::HasDof(DofType dofType) const
{
    auto it = m_dofInfo.find(dofType);
    return (it != m_dofInfo.end() && it->second.isActive);
}

//-----------------------------------------------------------------------------
int RgDofSchema::GetDofIndex(DofType dofType) const
{
    auto it = m_typeToIndex.find(dofType);
    return (it != m_typeToIndex.end()) ? it->second : -1;
}

//-----------------------------------------------------------------------------
int RgDofSchema::GetDofIndex(const std::string& varName) const
{
    auto it = m_nameToType.find(varName);
    if (it != m_nameToType.end())
    {
        return GetDofIndex(it->second);
    }
    return -1;
}

//-----------------------------------------------------------------------------
DofType RgDofSchema::GetDofType(int systemIndex) const
{
    if (systemIndex >= 0 && systemIndex < (int)m_activeDofs.size())
    {
        return m_activeDofs[systemIndex];
    }
    return DofType::MAX_DOF_TYPES;  // Invalid
}

//-----------------------------------------------------------------------------
const RgDofInfo* RgDofSchema::GetDofInfo(DofType dofType) const
{
    auto it = m_dofInfo.find(dofType);
    return (it != m_dofInfo.end()) ? &it->second : nullptr;
}

//-----------------------------------------------------------------------------
const RgDofInfo* RgDofSchema::GetDofInfo(int systemIndex) const
{
    if (systemIndex >= 0 && systemIndex < (int)m_activeDofs.size())
    {
        return GetDofInfo(m_activeDofs[systemIndex]);
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
std::vector<DofType> RgDofSchema::GetActiveDofs() const
{
    return m_activeDofs;
}

//-----------------------------------------------------------------------------
std::vector<RgDofInfo> RgDofSchema::GetAllDofInfo() const
{
    std::vector<RgDofInfo> result;
    for (DofType dof : m_activeDofs)
    {
        auto it = m_dofInfo.find(dof);
        if (it != m_dofInfo.end())
        {
            result.push_back(it->second);
        }
    }
    return result;
}

//-----------------------------------------------------------------------------
bool RgDofSchema::HasStructuralDofs() const
{
    return HasDof(DofType::DISPLACEMENT_X) ||
        HasDof(DofType::DISPLACEMENT_Y) ||
        HasDof(DofType::DISPLACEMENT_Z);
}

//-----------------------------------------------------------------------------
bool RgDofSchema::HasRotationalDofs() const
{
    return HasDof(DofType::ROTATION_X) ||
        HasDof(DofType::ROTATION_Y) ||
        HasDof(DofType::ROTATION_Z);
}

//-----------------------------------------------------------------------------
bool RgDofSchema::HasThermalDofs() const
{
    return HasDof(DofType::TEMPERATURE);
}

//-----------------------------------------------------------------------------
bool RgDofSchema::Is2D() const
{
    return HasDof(DofType::DISPLACEMENT_X) &&
        HasDof(DofType::DISPLACEMENT_Y) &&
        !HasDof(DofType::DISPLACEMENT_Z);
}

//-----------------------------------------------------------------------------
bool RgDofSchema::Is3D() const
{
    return HasDof(DofType::DISPLACEMENT_X) &&
        HasDof(DofType::DISPLACEMENT_Y) &&
        HasDof(DofType::DISPLACEMENT_Z);
}

//-----------------------------------------------------------------------------
int RgDofSchema::GetDimension() const
{
    if (Is3D()) return 3;
    if (Is2D()) return 2;
    return 0;  // Unclear or non-spatial
}

//-----------------------------------------------------------------------------
bool RgDofSchema::IsCompatible(const RgDofSet& dofSet) const
{
    // Check if all DOFs in the set are active
    const std::vector<DofType>& dofs = dofSet.GetDofs();
    for (DofType dof : dofs)
    {
        if (!HasDof(dof))
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
std::vector<int> RgDofSchema::GetIndices(const RgDofSet& dofSet) const
{
    std::vector<int> indices;
    const std::vector<DofType>& dofs = dofSet.GetDofs();

    for (DofType dof : dofs)
    {
        int idx = GetDofIndex(dof);
        if (idx >= 0)
        {
            indices.push_back(idx);
        }
    }

    return indices;
}

//-----------------------------------------------------------------------------
// Static helper functions
//-----------------------------------------------------------------------------

DofType RgDofSchema::GetDofTypeFromName(const std::string& name)
{
    if (name == DofVar::DISPLACEMENT_X) return DofType::DISPLACEMENT_X;
    if (name == DofVar::DISPLACEMENT_Y) return DofType::DISPLACEMENT_Y;
    if (name == DofVar::DISPLACEMENT_Z) return DofType::DISPLACEMENT_Z;
    if (name == DofVar::ROTATION_X) return DofType::ROTATION_X;
    if (name == DofVar::ROTATION_Y) return DofType::ROTATION_Y;
    if (name == DofVar::ROTATION_Z) return DofType::ROTATION_Z;
    if (name == DofVar::TEMPERATURE) return DofType::TEMPERATURE;
    if (name == DofVar::VELOCITY_X) return DofType::VELOCITY_X;
    if (name == DofVar::VELOCITY_Y) return DofType::VELOCITY_Y;
    if (name == DofVar::VELOCITY_Z) return DofType::VELOCITY_Z;
    if (name == DofVar::PRESSURE) return DofType::PRESSURE;

    return DofType::MAX_DOF_TYPES;  // Invalid
}

//-----------------------------------------------------------------------------
std::string RgDofSchema::GetDofName(DofType dofType)
{
    switch (dofType)
    {
    case DofType::DISPLACEMENT_X: return DofVar::DISPLACEMENT_X;
    case DofType::DISPLACEMENT_Y: return DofVar::DISPLACEMENT_Y;
    case DofType::DISPLACEMENT_Z: return DofVar::DISPLACEMENT_Z;
    case DofType::ROTATION_X: return DofVar::ROTATION_X;
    case DofType::ROTATION_Y: return DofVar::ROTATION_Y;
    case DofType::ROTATION_Z: return DofVar::ROTATION_Z;
    case DofType::TEMPERATURE: return DofVar::TEMPERATURE;
    case DofType::VELOCITY_X: return DofVar::VELOCITY_X;
    case DofType::VELOCITY_Y: return DofVar::VELOCITY_Y;
    case DofType::VELOCITY_Z: return DofVar::VELOCITY_Z;
    case DofType::PRESSURE: return DofVar::PRESSURE;
    default: return "unknown";
    }
}

//-----------------------------------------------------------------------------
std::string RgDofSchema::GetDofDescription(DofType dofType)
{
    switch (dofType)
    {
    case DofType::DISPLACEMENT_X: return "X displacement";
    case DofType::DISPLACEMENT_Y: return "Y displacement";
    case DofType::DISPLACEMENT_Z: return "Z displacement";
    case DofType::ROTATION_X: return "Rotation about X";
    case DofType::ROTATION_Y: return "Rotation about Y";
    case DofType::ROTATION_Z: return "Rotation about Z";
    case DofType::TEMPERATURE: return "Temperature";
    case DofType::VELOCITY_X: return "X velocity";
    case DofType::VELOCITY_Y: return "Y velocity";
    case DofType::VELOCITY_Z: return "Z velocity";
    case DofType::PRESSURE: return "Pressure";
    default: return "Unknown DOF";
    }
}

//-----------------------------------------------------------------------------
void RgDofSchema::PrintSchema() const
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("DOF Schema Configuration\n");
    RgLog("=================================================\n");
    RgLog("DOFs per node: %d\n", GetDofsPerNode());
    RgLog("-------------------------------------------------\n");
    RgLog("Index  Name    Description\n");
    RgLog("-------------------------------------------------\n");

    for (int i = 0; i < (int)m_activeDofs.size(); ++i)
    {
        const RgDofInfo* info = GetDofInfo(i);
        if (info)
        {
            RgLog("  %2d   %-6s  %s\n",
                info->systemIndex,
                info->name.c_str(),
                info->description.c_str());
        }
    }

    RgLog("=================================================\n\n");
}

//-----------------------------------------------------------------------------
bool RgDofSchema::Validate() const
{
    // Check for duplicate indices
    std::set<int> indices;
    for (const auto& dof : m_activeDofs)
    {
        int idx = GetDofIndex(dof);
        if (indices.find(idx) != indices.end())
        {
            RgLogError("Duplicate DOF index: %d", idx);
            return false;
        }
        indices.insert(idx);
    }

    // Check index continuity (should be 0, 1, 2, ...)
    for (int i = 0; i < (int)m_activeDofs.size(); ++i)
    {
        int idx = GetDofIndex(m_activeDofs[i]);
        if (idx != i)
        {
            RgLogError("DOF index not continuous: expected %d, got %d", i, idx);
            return false;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
void RgDofSchema::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        // Save number of active DOFs
        int nDofs = GetDofsPerNode();
        ar << nDofs;

        // Save each active DOF type
        for (DofType dof : m_activeDofs)
        {
            int dofInt = static_cast<int>(dof);
            ar << dofInt;
        }
    }
    else
    {
        Clear();

        // Load number of DOFs
        int nDofs;
        ar >> nDofs;

        // Load each DOF type and activate it
        for (int i = 0; i < nDofs; ++i)
        {
            int dofInt;
            ar >> dofInt;
            DofType dof = static_cast<DofType>(dofInt);
            AddDof(dof);
        }
    }
}

//=============================================================================
// RgDofSchemaPresets Implementation
//=============================================================================

RgDofSchema* RgDofSchemaPresets::CreateStructural3D()
{
    RgDofSchema* schema = new RgDofSchema();
    schema->SetSchema(RgDofSet::Solid3D());
    return schema;
}

RgDofSchema* RgDofSchemaPresets::CreateStructural2D()
{
    RgDofSchema* schema = new RgDofSchema();
    schema->SetSchema(RgDofSet::Plane2D());
    return schema;
}

RgDofSchema* RgDofSchemaPresets::CreateBeamShell3D()
{
    RgDofSchema* schema = new RgDofSchema();
    schema->SetSchema(RgDofSet::Beam3D());
    return schema;
}

RgDofSchema* RgDofSchemaPresets::CreateThermal()
{
    RgDofSchema* schema = new RgDofSchema();
    schema->SetSchema(RgDofSet::Thermal3D());
    return schema;
}

RgDofSchema* RgDofSchemaPresets::CreateThermoMechanical()
{
    RgDofSchema* schema = new RgDofSchema();
    schema->SetSchema(RgDofSet::Coupled_ThermoMech());
    return schema;
}

RgDofSchema* RgDofSchemaPresets::CreateFluidDynamics()
{
    RgDofSchema* schema = new RgDofSchema();
    RgDofSet fluidSet;
    fluidSet.AddDof(DofType::VELOCITY_X);
    fluidSet.AddDof(DofType::VELOCITY_Y);
    fluidSet.AddDof(DofType::VELOCITY_Z);
    fluidSet.AddDof(DofType::PRESSURE);
    schema->SetSchema(fluidSet);
    return schema;
}