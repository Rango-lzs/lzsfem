/*********************************************************************
 * \file   RgDofSchema.h
 * \brief  DOF Schema - Defines which DOFs are active in the system
 *
 * This class describes the CONFIGURATION of DOFs in the FE system,
 * not the management of individual node DOFs.
 *
 * Naming Convention:
 * - RgDofSchema: Describes the DOF configuration/blueprint
 * - RgDofInfo: Information about a single DOF type
 * - RgDofSet: A collection of DOF types (for element types)
 * - FENode: Actual DOF values at a node
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#pragma once

#include "femcore/fem_export.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

 //=============================================================================
 // DOF Type Enumerations
 //=============================================================================

 /// Basic DOF types for different physics
enum class DofType
{
    // Structural DOFs (0-5)
    DISPLACEMENT_X = 0,  ///< u - displacement in X
    DISPLACEMENT_Y = 1,  ///< v - displacement in Y
    DISPLACEMENT_Z = 2,  ///< w - displacement in Z
    ROTATION_X = 3,      ///< θx - rotation about X
    ROTATION_Y = 4,      ///< θy - rotation about Y
    ROTATION_Z = 5,      ///< θz - rotation about Z

    // Thermal DOFs (6-6)
    TEMPERATURE = 6,     ///< T - temperature

    // Fluid DOFs (7-10)
    VELOCITY_X = 7,      ///< vx - fluid velocity in X
    VELOCITY_Y = 8,      ///< vy - fluid velocity in Y
    VELOCITY_Z = 9,      ///< vz - fluid velocity in Z
    PRESSURE = 10,       ///< p - fluid pressure

    // Electromagnetic DOFs (11-13)
    ELECTRIC_POTENTIAL = 11,  ///< φ - electric potential
    MAGNETIC_POTENTIAL_X = 12, ///< Ax - magnetic vector potential X
    MAGNETIC_POTENTIAL_Y = 13, ///< Ay - magnetic vector potential Y
    MAGNETIC_POTENTIAL_Z = 14, ///< Az - magnetic vector potential Z

    // Reserved for future use (15-31)
    CUSTOM_1 = 15,
    CUSTOM_2 = 16,
    CUSTOM_3 = 17,
    CUSTOM_4 = 18,

    MAX_DOF_TYPES = 32   ///< Maximum number of DOF types
};

/// DOF variable names (for string-based access)
namespace DofVar
{
    // Structural
    constexpr const char* DISPLACEMENT_X = "u";
    constexpr const char* DISPLACEMENT_Y = "v";
    constexpr const char* DISPLACEMENT_Z = "w";
    constexpr const char* ROTATION_X = "Rx";
    constexpr const char* ROTATION_Y = "Ry";
    constexpr const char* ROTATION_Z = "Rz";

    // Thermal
    constexpr const char* TEMPERATURE = "T";

    // Fluid
    constexpr const char* VELOCITY_X = "vx";
    constexpr const char* VELOCITY_Y = "vy";
    constexpr const char* VELOCITY_Z = "vz";
    constexpr const char* PRESSURE = "p";

    // Electromagnetic
    constexpr const char* ELECTRIC_POTENTIAL = "phi";
    constexpr const char* MAGNETIC_POTENTIAL_X = "Ax";
    constexpr const char* MAGNETIC_POTENTIAL_Y = "Ay";
    constexpr const char* MAGNETIC_POTENTIAL_Z = "Az";
}

//=============================================================================
// RgDofInfo - Information about a single DOF type
//=============================================================================

/// Information about a single degree of freedom
struct FEM_EXPORT RgDofInfo
{
    DofType type;           ///< DOF type enumeration
    std::string name;       ///< Variable name (e.g., "u", "v", "T")
    std::string description; ///< Full description
    int systemIndex;        ///< Index in the system (assigned by schema)
    bool isActive;          ///< Whether this DOF is active in current analysis

    RgDofInfo()
        : type(DofType::DISPLACEMENT_X)
        , name("")
        , description("")
        , systemIndex(-1)
        , isActive(false)
    {}

    RgDofInfo(DofType t, const std::string& n, const std::string& d = "")
        : type(t)
        , name(n)
        , description(d)
        , systemIndex(-1)
        , isActive(false)
    {}
};

//=============================================================================
// RgDofSet - A collection of DOFs (e.g., for a specific element type)
//=============================================================================

/// Represents a set of DOFs (for example, what a specific element type uses)
class FEM_EXPORT RgDofSet
{
public:
    RgDofSet() {}

    /// Add a DOF to this set
    void AddDof(DofType dof);

    /// Check if DOF is in this set
    bool HasDof(DofType dof) const;

    /// Get number of DOFs in this set
    int Size() const { return (int)m_dofs.size(); }

    /// Get DOF at index
    DofType GetDof(int index) const;

    /// Get all DOFs
    const std::vector<DofType>& GetDofs() const { return m_dofs; }

    /// Clear all DOFs
    void Clear() { m_dofs.clear(); }

    void Sort() { std::sort(m_dofs.begin(), m_dofs.end()); }

    /// Create predefined DOF sets
    static RgDofSet Solid3D();      ///< 3D solid: u, v, w
    static RgDofSet Shell();        ///< Shell: u, v, w, Rx, Ry, Rz
    static RgDofSet Beam3D();       ///< 3D beam: u, v, w, Rx, Ry, Rz
    static RgDofSet Truss3D();      ///< 3D truss: u, v, w
    static RgDofSet Plane2D();      ///< 2D plane: u, v
    static RgDofSet Thermal3D();    ///< 3D thermal: T
    static RgDofSet Coupled_ThermoMech(); ///< Thermo-mechanical: u,v,w,T

private:
    std::vector<DofType> m_dofs;
};

//=============================================================================
// RgDofSchema - System DOF Configuration/Blueprint
//=============================================================================

/// Describes which DOFs are active in the FE system
/// This is the "schema" or "configuration" that defines the structure
/// of the equation system, not the manager of node DOF values.
class FEM_EXPORT RgDofSchema
{
public:
    RgDofSchema();
    ~RgDofSchema();

    //-------------------------------------------------------------------------
    // Schema Setup
    //-------------------------------------------------------------------------

    /// Initialize schema with a predefined DOF set
    void SetSchema(const RgDofSet& dofSet);

    /// Add a single DOF type to the schema
    bool AddDof(DofType dofType, const std::string& name = "",
        const std::string& description = "");

    /// Add multiple DOF types at once
    bool AddDofs(const std::vector<DofType>& dofTypes);

    /// Remove a DOF type from the schema
    bool RemoveDof(DofType dofType);

    /// Clear all DOFs
    void Clear();

    /// Reset the schema (same as Clear)
    void Reset() { Clear(); }

    //-------------------------------------------------------------------------
    // Schema Query
    //-------------------------------------------------------------------------

    /// Get total number of DOFs per node in this schema
    int GetDofsPerNode() const { return (int)m_activeDofs.size(); }

    /// Check if a DOF type is in this schema
    bool HasDof(DofType dofType) const;

    /// Get system index for a DOF type (-1 if not in schema)
    int GetDofIndex(DofType dofType) const;

    /// Get system index by variable name (e.g., "u", "v", "T")
    int GetDofIndex(const std::string& varName) const;

    /// Get DOF type from system index
    DofType GetDofType(int systemIndex) const;

    /// Get DOF info for a specific type
    const RgDofInfo* GetDofInfo(DofType dofType) const;

    /// Get DOF info by system index
    const RgDofInfo* GetDofInfo(int systemIndex) const;

    /// Get all active DOF types
    std::vector<DofType> GetActiveDofs() const;

    /// Get all DOF information
    std::vector<RgDofInfo> GetAllDofInfo() const;

    //-------------------------------------------------------------------------
    // Physics-Specific Queries
    //-------------------------------------------------------------------------

    /// Check if structural DOFs are in schema (u,v,w)
    bool HasStructuralDofs() const;

    /// Check if rotational DOFs are in schema (Rx,Ry,Rz)
    bool HasRotationalDofs() const;

    /// Check if thermal DOFs are in schema (T)
    bool HasThermalDofs() const;

    /// Check if this is a 2D schema (only u,v active)
    bool Is2D() const;

    /// Check if this is a 3D schema (u,v,w active)
    bool Is3D() const;

    /// Get spatial dimensionality (2 or 3)
    int GetDimension() const;

    //-------------------------------------------------------------------------
    // Element Compatibility
    //-------------------------------------------------------------------------

    /// Check if a DOF set is compatible with this schema
    /// (all DOFs in the set are present in the schema)
    bool IsCompatible(const RgDofSet& dofSet) const;

    /// Get system indices for DOFs in a set
    std::vector<int> GetIndices(const RgDofSet& dofSet) const;

    //-------------------------------------------------------------------------
    // String-based Access (for input files, scripts, etc.)
    //-------------------------------------------------------------------------

    /// Get DOF type from variable name
    static DofType GetDofTypeFromName(const std::string& name);

    /// Get variable name from DOF type
    static std::string GetDofName(DofType dofType);

    /// Get description from DOF type
    static std::string GetDofDescription(DofType dofType);

    //-------------------------------------------------------------------------
    // Utility Functions
    //-------------------------------------------------------------------------

    /// Print schema information to log
    void PrintSchema() const;

    /// Validate schema configuration
    bool Validate() const;

    /// Serialize to/from restart file
    void Serialize(class DumpStream& ar);

private:
    /// Map of all DOF information (indexed by DofType)
    std::map<DofType, RgDofInfo> m_dofInfo;

    /// List of active DOFs (in order of system index)
    std::vector<DofType> m_activeDofs;

    /// Map from variable name to DOF type
    std::map<std::string, DofType> m_nameToType;

    /// Map from DOF type to system index
    std::map<DofType, int> m_typeToIndex;

    /// Rebuild index maps after changes
    void RebuildIndexMaps();

    /// Create default DOF information
    void CreateDefaultDofInfo();
};

//=============================================================================
// Helper Functions - Schema Presets
//=============================================================================

/// Create DOF schemas for common analysis types
namespace RgDofSchemaPresets
{
    /// 3D structural analysis (u,v,w)
    FEM_EXPORT RgDofSchema* CreateStructural3D();

    /// 2D structural analysis (u,v)
    FEM_EXPORT RgDofSchema* CreateStructural2D();

    /// 3D beam/shell analysis (u,v,w,Rx,Ry,Rz)
    FEM_EXPORT RgDofSchema* CreateBeamShell3D();

    /// Thermal analysis (T)
    FEM_EXPORT RgDofSchema* CreateThermal();

    /// Coupled thermo-mechanical (u,v,w,T)
    FEM_EXPORT RgDofSchema* CreateThermoMechanical();

    /// Fluid dynamics (vx,vy,vz,p)
    FEM_EXPORT RgDofSchema* CreateFluidDynamics();
}