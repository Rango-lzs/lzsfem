// ============================================================================
// FEM Data Kernel - Abaqus-like Design
// ============================================================================

#ifndef FEM_KERNEL_H
#define FEM_KERNEL_H

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
// Node Class
// ============================================================================

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

// ============================================================================
// Element Class
// ============================================================================

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

// ============================================================================
// Sets for Grouping
// ============================================================================

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

// ============================================================================
// Material Class
// ============================================================================

class Material {
public:
    explicit Material(const std::string& name) : name_(name) {}
    
    void setElastic(Real E, Real nu) {
        properties_["E"] = E;
        properties_["nu"] = nu;
    }
    
    void setDensity(Real rho) {
        properties_["density"] = rho;
    }
    
    void setProperty(const std::string& key, Real value) {
        properties_[key] = value;
    }
    
    Real getProperty(const std::string& key) const {
        auto it = properties_.find(key);
        return (it != properties_.end()) ? it->second : 0.0;
    }
    
    const std::string& getName() const { return name_; }

private:
    std::string name_;
    std::map<std::string, Real> properties_;
};

// ============================================================================
// Section Class
// ============================================================================

enum class SectionType {
    SOLID,
    SHELL,
    BEAM
};

class Section {
public:
    Section(const std::string& name, SectionType type)
        : name_(name), type_(type) {}
    
    void setMaterial(int matId) { materialId_ = matId; }
    int getMaterial() const { return materialId_; }
    
    void setThickness(Real t) { thickness_ = t; }  // For shell
    Real getThickness() const { return thickness_; }
    
    void setProperty(const std::string& key, Real value) {
        properties_[key] = value;
    }
    
    const std::string& getName() const { return name_; }
    SectionType getType() const { return type_; }

private:
    std::string name_;
    SectionType type_;
    int materialId_ = -1;
    Real thickness_ = 0.0;
    std::map<std::string, Real> properties_;
};

// ============================================================================
// Boundary Condition
// ============================================================================

enum class BCType {
    DISPLACEMENT,
    VELOCITY,
    ACCELERATION,
    TEMPERATURE
};

class BoundaryCondition {
public:
    BoundaryCondition(const std::string& name, BCType type)
        : name_(name), type_(type) {}
    
    void apply(const std::string& nodeSetName, size_t dof, Real value) {
        BCData data;
        data.nodeSetName = nodeSetName;
        data.dof = dof;
        data.value = value;
        conditions_.push_back(data);
    }
    
    struct BCData {
        std::string nodeSetName;
        size_t dof;
        Real value;
    };
    
    const std::vector<BCData>& getConditions() const { return conditions_; }
    const std::string& getName() const { return name_; }

private:
    std::string name_;
    BCType type_;
    std::vector<BCData> conditions_;
};

// ============================================================================
// Load Class
// ============================================================================

enum class LoadType {
    CONCENTRATED_FORCE,
    DISTRIBUTED_LOAD,
    PRESSURE,
    BODY_FORCE,
    THERMAL
};

class Load {
public:
    Load(const std::string& name, LoadType type)
        : name_(name), type_(type) {}
    
    void applyConcentratedForce(const std::string& nodeSetName, 
                                size_t dof, Real magnitude) {
        LoadData data;
        data.targetSetName = nodeSetName;
        data.dof = dof;
        data.magnitude = magnitude;
        loads_.push_back(data);
    }
    
    void applyPressure(const std::string& surfaceName, Real magnitude) {
        LoadData data;
        data.targetSetName = surfaceName;
        data.magnitude = magnitude;
        loads_.push_back(data);
    }
    
    struct LoadData {
        std::string targetSetName;
        size_t dof = 0;
        Real magnitude = 0.0;
        Vec3 direction = {0, 0, 0};
    };
    
    const std::vector<LoadData>& getLoads() const { return loads_; }
    const std::string& getName() const { return name_; }

private:
    std::string name_;
    LoadType type_;
    std::vector<LoadData> loads_;
};

// ============================================================================
// Part Class
// ============================================================================

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

// ============================================================================
// Model Class - Top Level Container
// ============================================================================

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
    void assembleModel() {
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

#endif // FEM_KERNEL_H