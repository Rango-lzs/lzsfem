// ============================================================================
// FEM Data Kernel - Model Class
// ============================================================================

#ifndef FEM_MODEL_H
#define FEM_MODEL_H

#include "BaseTypes.h"
#include "Part.h"
#include "Material.h"
#include "Section.h"
#include "BoundaryCondition.h"
#include "Load.h"

#include <map>
#include <memory>
#include <string>

namespace FEM {

class Model {
public:
    explicit Model(const std::string& name) : name_(name) {}
    
    // Part management
    void addPart(const Part& part) {
        parts_[part.getName()] = std::make_shared<Part>(part);
    }
    
    std::shared_ptr<Part> getPart(const std::string& name) {
        auto it = parts_.find(name);
        return (it != parts_.end()) ? it->second : nullptr;
    }
    
    // Material management
    int addMaterial(const Material& mat) {
        int id = static_cast<int>(materials_.size());
        materials_.push_back(std::make_shared<Material>(mat));
        materialMap_[mat.getName()] = id;
        return id;
    }
    
    std::shared_ptr<Material> getMaterial(int id) {
        return (id >= 0 && id < materials_.size()) ? materials_[id] : nullptr;
    }
    
    std::shared_ptr<Material> getMaterial(const std::string& name) {
        auto it = materialMap_.find(name);
        if (it != materialMap_.end()) {
            return materials_[it->second];
        }
        return nullptr;
    }
    
    // Section management
    int addSection(const Section& sec) {
        int id = static_cast<int>(sections_.size());
        sections_.push_back(std::make_shared<Section>(sec));
        sectionMap_[sec.getName()] = id;
        return id;
    }
    
    std::shared_ptr<Section> getSection(int id) {
        return (id >= 0 && id < sections_.size()) ? sections_[id] : nullptr;
    }
    
    // Boundary conditions and loads
    void addBoundaryCondition(const BoundaryCondition& bc) {
        boundaryConditions_[bc.getName()] = 
            std::make_shared<BoundaryCondition>(bc);
    }
    
    void addLoad(const Load& load) {
        loads_[load.getName()] = std::make_shared<Load>(load);
    }
    
    // Assembly - flattened mesh for analysis
    void assembleModel();
    
    const std::map<NodeID, std::shared_ptr<Node>>& getGlobalNodes() const {
        return globalNodes_;
    }
    
    const std::map<ElementID, std::shared_ptr<Element>>& getGlobalElements() const {
        return globalElements_;
    }
    
    const std::string& getName() const { return name_; }

private:
    std::string name_;
    
    // Parts
    std::map<std::string, std::shared_ptr<Part>> parts_;
    
    // Materials and sections
    std::vector<std::shared_ptr<Material>> materials_;
    std::map<std::string, int> materialMap_;
    std::vector<std::shared_ptr<Section>> sections_;
    std::map<std::string, int> sectionMap_;
    
    // Boundary conditions and loads
    std::map<std::string, std::shared_ptr<BoundaryCondition>> boundaryConditions_;
    std::map<std::string, std::shared_ptr<Load>> loads_;
    
    // Assembled global mesh
    std::map<NodeID, std::shared_ptr<Node>> globalNodes_;
    std::map<ElementID, std::shared_ptr<Element>> globalElements_;
};

} // namespace FEM

#endif // FEM_MODEL_H