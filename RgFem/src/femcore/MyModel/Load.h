// ============================================================================
// FEM Data Kernel - Load Class
// ============================================================================

#ifndef FEM_LOAD_H
#define FEM_LOAD_H

#include "BaseTypes.h"
#include <string>
#include <vector>

namespace FEM {

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

} // namespace FEM

#endif // FEM_LOAD_H