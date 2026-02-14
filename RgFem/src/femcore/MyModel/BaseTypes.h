// ============================================================================
// FEM Data Kernel - Base Types
// ============================================================================

#ifndef FEM_BASE_TYPES_H
#define FEM_BASE_TYPES_H

#include <vector>
#include <map>
#include <string>
#include <memory>
#include <set>
#include <array>

namespace FEM {

// ============================================================================
// Basic Types
// ============================================================================

using NodeID = int;
using ElementID = int;
using Real = double;

template<size_t N>
using Vector = std::array<Real, N>;

using Vec3 = Vector<3>;

// ============================================================================
// Element Types
// ============================================================================

enum class ElementType {
    C3D4,   // 4-node tetrahedron
    C3D8,   // 8-node brick
    C3D10,  // 10-node tetrahedron
    C3D20,  // 20-node brick
    S3,     // 3-node triangle shell
    S4,     // 4-node quad shell
    B31,    // 2-node beam
    SPRING, // Spring element
    CUSTOM
};

enum class SectionType {
    SOLID,
    SHELL,
    BEAM
};

enum class BCType {
    DISPLACEMENT,
    VELOCITY,
    ACCELERATION,
    TEMPERATURE
};

enum class LoadType {
    CONCENTRATED_FORCE,
    DISTRIBUTED_LOAD,
    PRESSURE,
    BODY_FORCE,
    THERMAL
};


} // namespace FEM

#endif // FEM_BASE_TYPES_H