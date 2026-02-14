// ============================================================================
// FEM Data Kernel - Section Class
// ============================================================================

#ifndef FEM_SECTION_H
#define FEM_SECTION_H

#include "BaseTypes.h"
#include <map>
#include <string>

namespace FEM {

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

} // namespace FEM

#endif // FEM_SECTION_H