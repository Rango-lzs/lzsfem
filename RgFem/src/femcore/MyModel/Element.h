// ============================================================================
// FEM Data Kernel - Element Class
// ============================================================================

#ifndef FEM_ELEMENT_H
#define FEM_ELEMENT_H

#include "BaseTypes.h"
#include <vector>

namespace FEM {

class Element {
public:
    Element(ElementID id, ElementType type, const std::vector<NodeID>& connectivity)
        : id_(id), type_(type), connectivity_(connectivity) {}
    
    ElementID getId() const { return id_; }
    ElementType getType() const { return type_; }
    const std::vector<NodeID>& getConnectivity() const { return connectivity_; }
    size_t getNumNodes() const { return connectivity_.size(); }
    
    // Property assignment
    void setMaterialId(int matId) { materialId_ = matId; }
    void setSectionId(int secId) { sectionId_ = secId; }
    int getMaterialId() const { return materialId_; }
    int getSectionId() const { return sectionId_; }

private:
    ElementID id_;
    ElementType type_;
    std::vector<NodeID> connectivity_;
    int materialId_ = -1;
    int sectionId_ = -1;
};

} // namespace FEM

#endif // FEM_ELEMENT_H