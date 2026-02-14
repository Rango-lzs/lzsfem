// ============================================================================
// FEM Data Kernel - Node Class
// ============================================================================

#ifndef FEM_NODE_H
#define FEM_NODE_H

#include "BaseTypes.h"

namespace FEM {

class Node {
public:
    Node(NodeID id, const Vec3& coords) 
        : id_(id), coords_(coords), dofs_(3, 0.0) {}
    
    NodeID getId() const { return id_; }
    const Vec3& getCoords() const { return coords_; }
    void setCoords(const Vec3& coords) { coords_ = coords; }
    
    // DOF management
    size_t getNumDofs() const { return dofs_.size(); }
    void setNumDofs(size_t n) { dofs_.resize(n, 0.0); }
    Real getDof(size_t i) const { return dofs_[i]; }
    void setDof(size_t i, Real value) { dofs_[i] = value; }

private:
    NodeID id_;
    Vec3 coords_;
    std::vector<Real> dofs_;  // Displacement, rotation, temperature, etc.
};

} // namespace FEM

#endif // FEM_NODE_H