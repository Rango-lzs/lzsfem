// ============================================================================
// FEM Data Kernel - Element Set Class
// ============================================================================

#ifndef FEM_ELEMENTSET_H
#define FEM_ELEMENTSET_H

#include "BaseTypes.h"
#include <set>
#include <string>
#include <vector>

namespace FEM {

class ElementSet {
public:
    explicit ElementSet(const std::string& name) : name_(name) {}
    
    void addElement(ElementID id) { elements_.insert(id); }
    void addElements(const std::vector<ElementID>& ids) {
        elements_.insert(ids.begin(), ids.end());
    }
    
    const std::set<ElementID>& getElements() const { return elements_; }
    const std::string& getName() const { return name_; }
    size_t size() const { return elements_.size(); }

private:
    std::string name_;
    std::set<ElementID> elements_;
};

} // namespace FEM

#endif // FEM_ELEMENTSET_H