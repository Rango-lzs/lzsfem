// ============================================================================
// FEM Data Kernel - Material Class
// ============================================================================

#ifndef FEM_MATERIAL_H
#define FEM_MATERIAL_H

#include "BaseTypes.h"
#include <map>
#include <string>

namespace FEM {

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

} // namespace FEM

#endif // FEM_MATERIAL_H