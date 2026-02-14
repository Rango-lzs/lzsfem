// ============================================================================
// FEM Data Kernel - Part Class
// ============================================================================

#ifndef FEM_PART_H
#define FEM_PART_H

#include "BaseTypes.h"
#include "Node.h"
#include "Element.h"
#include "NodeSet.h"
#include "ElementSet.h"

#include <map>
#include <memory>
#include <string>

namespace FEM {

class Part {
public:
    explicit Part(const std::string& name) : name_(name) {}
    
    // Node management
    void addNode(const Node& node) {
        nodes_[node.getId()] = std::make_shared<Node>(node);
    }
    
    std::shared_ptr<Node> getNode(NodeID id) {
        auto it = nodes_.find(id);
        return (it != nodes_.end()) ? it->second : nullptr;
    }
    
    const std::map<NodeID, std::shared_ptr<Node>>& getNodes() const {
        return nodes_;
    }
    
    // Element management
    void addElement(const Element& elem) {
        elements_[elem.getId()] = std::make_shared<Element>(elem);
    }
    
    std::shared_ptr<Element> getElement(ElementID id) {
        auto it = elements_.find(id);
        return (it != elements_.end()) ? it->second : nullptr;
    }
    
    const std::map<ElementID, std::shared_ptr<Element>>& getElements() const {
        return elements_;
    }
    
    // Set management
    void addNodeSet(const NodeSet& nset) {
        nodeSets_[nset.getName()] = std::make_shared<NodeSet>(nset);
    }
    
    void addElementSet(const ElementSet& eset) {
        elementSets_[eset.getName()] = std::make_shared<ElementSet>(eset);
    }
    
    std::shared_ptr<NodeSet> getNodeSet(const std::string& name) {
        auto it = nodeSets_.find(name);
        return (it != nodeSets_.end()) ? it->second : nullptr;
    }
    
    std::shared_ptr<ElementSet> getElementSet(const std::string& name) {
        auto it = elementSets_.find(name);
        return (it != elementSets_.end()) ? it->second : nullptr;
    }
    
    const std::string& getName() const { return name_; }

private:
    std::string name_;
    std::map<NodeID, std::shared_ptr<Node>> nodes_;
    std::map<ElementID, std::shared_ptr<Element>> elements_;
    std::map<std::string, std::shared_ptr<NodeSet>> nodeSets_;
    std::map<std::string, std::shared_ptr<ElementSet>> elementSets_;
};

} // namespace FEM

#endif // FEM_PART_H