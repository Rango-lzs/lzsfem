// ============================================================================
// FEM Data Kernel - Node Set Class
// ============================================================================

#ifndef FEM_NODESET_H
#define FEM_NODESET_H

#include "BaseTypes.h"
#include <set>
#include <string>

namespace FEM {

class NodeSet {
public:
    explicit NodeSet(const std::string& name) : name_(name) {}
    
    void addNode(NodeID id) { nodes_.insert(id); }
    void addNodes(const std::vector<NodeID>& ids) {
        nodes_.insert(ids.begin(), ids.end());
    }
    
    const std::set<NodeID>& getNodes() const { return nodes_; }
    const std::string& getName() const { return name_; }
    size_t size() const { return nodes_.size(); }

private:
    std::string name_;
    std::set<NodeID> nodes_;
};

} // namespace FEM

#endif // FEM_NODESET_H