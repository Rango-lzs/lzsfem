// ============================================================================
// FEM Data Kernel - Model Class Implementation
// ============================================================================

#include "Model.h"
#include "Node.h"

namespace FEM {

void Model::assembleModel() {
    // Flatten all parts into global nodes and elements
    NodeID globalNodeId = 1;
    ElementID globalElemId = 1;
    
    for (auto& [partName, part] : parts_) {
        // Copy nodes with new IDs
        for (auto& [localId, node] : part->getNodes()) {
            Node globalNode(globalNodeId++, node->getCoords());
            globalNodes_[globalNode.getId()] = 
                std::make_shared<Node>(globalNode);
        }
        
        // Copy elements with new IDs and updated connectivity
        for (auto& [localId, elem] : part->getElements()) {
            Element globalElem(globalElemId++, elem->getType(), 
                               elem->getConnectivity());
            globalElem.setMaterialId(elem->getMaterialId());
            globalElem.setSectionId(elem->getSectionId());
            globalElements_[globalElem.getId()] = 
                std::make_shared<Element>(globalElem);
        }
    }
}

} // namespace FEM