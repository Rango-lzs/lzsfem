// ============================================================================
// FEM Data Kernel - Boundary Condition Class
// ============================================================================

#ifndef FEM_BOUNDARYCONDITION_H
#define FEM_BOUNDARYCONDITION_H

#include "BaseTypes.h"
#include <string>
#include <vector>

namespace FEM {

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

} // namespace FEM

#endif // FEM_BOUNDARYCONDITION_H