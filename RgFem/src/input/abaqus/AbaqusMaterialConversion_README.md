/*********************************************************************
 * \file   AbaqusMaterialConversion_README.md
 * \brief  Guide for Abaqus material conversion implementation
 *
 * \author
 * \date   February 2026
 *********************************************************************/

# Abaqus Material Conversion Implementation Guide

## Overview

This guide explains how to integrate the Abaqus material conversion functionality
into your FE solver. The implementation converts Abaqus material definitions from
INP files into your FEModel material objects.

## Files Added

1. **AbaqusMaterialConverter.h** - Header for material conversion utilities
2. **AbaqusMaterialConverter.cpp** - Implementation of conversion logic
3. **AbaqusImport_MaterialUpdate.cpp** - Updated parsing and processing functions

## Integration Steps

### Step 1: Add Headers to AbaqusImport

Add these includes to `AbaqusImport.cpp`:

```cpp
#include "AbaqusMaterialConverter.h"
#include "materials/SmallDeformation/RgLinearElasticMaterial.h"
#include "materials/SmallDeformation/RgElastoPlasticMaterial.h"
```

### Step 2: Update AbaqusImport.h

Add the material creation method to the `AbaqusImport` class:

```cpp
private:
    // ... existing members ...
    
    /**
     * @brief Create FEModel materials from parsed Abaqus materials
     */
    bool createMaterials(FEModel* fem);
```

### Step 3: Replace parseMaterial Function

Replace the existing `parseMaterial` function in `AbaqusImport.cpp` with the
updated version from `AbaqusImport_MaterialUpdate.cpp`. This new version:

- Handles ELASTIC properties (E, nu)
- Handles PLASTIC properties (yield stress, hardening curve)
- Handles DENSITY
- Provides better error handling and logging

### Step 4: Add createMaterials Function

Add the `createMaterials` function from `AbaqusImport_MaterialUpdate.cpp` to
`AbaqusImport.cpp`. This function:

- Iterates through all parsed materials
- Uses `AbaqusMaterialConverter` to create FEModel materials
- Adds materials to the FEModel
- Provides comprehensive logging

### Step 5: Update processData Function

In the `processData` function, add material creation BEFORE processing parts:

```cpp
bool AbaqusImport::processData(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();
    
    // 1. Create materials FIRST
    if (!createMaterials(fem))
    {
        RgLogError("Failed to create materials");
        return false;
    }
    
    // 2. Then process geometry (existing code)
    // ...
}
```

### Step 6: Optional - Add Material Base Class Methods

If your `RgMaterial` class doesn't have these methods, consider adding them:

```cpp
class RgMaterial
{
public:
    // ... existing methods ...
    
    // Optional: for material identification
    void SetID(int id) { m_id = id; }
    int GetID() const { return m_id; }
    
    // Optional: for material naming
    void SetName(const std::string& name) { m_name = name; }
    std::string GetName() const { return m_name; }
    
    // For density handling
    void SetDensity(double rho) { m_density = rho; }
    double GetDensity() const { return m_density; }
    
private:
    int m_id;
    std::string m_name;
    double m_density;
};
```

## Supported Material Types

### 1. Linear Elastic

**Abaqus INP format:**
```
*Material, name=Steel
*Elastic
210000, 0.3
*Density
7850
```

**Converts to:** `SmallDef::RgLinearElastic(E=210000, nu=0.3)`

### 2. Elasto-Plastic (Linear Hardening)

**Abaqus INP format:**
```
*Material, name=Aluminum
*Elastic
70000, 0.33
*Plastic
250, 0.0
350, 0.1
*Density
2700
```

**Converts to:** `SmallDef::RgElastoPlastic(E=70000, nu=0.33, sy=250, H=1000)`

The hardening modulus H is calculated from the plastic curve:
- H = (stress2 - yield_stress) / plastic_strain2
- For the example: H = (350 - 250) / 0.1 = 1000

### 3. Perfectly Plastic

**Abaqus INP format:**
```
*Material, name=SoftSteel
*Elastic
200000, 0.3
*Plastic
200
```

**Converts to:** `SmallDef::RgElastoPlastic(E=200000, nu=0.3, sy=200, H=0)`

## Material-Domain Assignment

### Current Implementation

Currently, materials are created but not automatically assigned to domains.
You need to implement section parsing to establish the material-domain mapping.

### Recommended Implementation

Add section data structure to `AbaqusImport`:

```cpp
struct SectionInfo
{
    std::string elementSetName;  // Which elements
    std::string materialName;     // Which material
    std::string type;             // SOLID, SHELL, etc.
};

std::vector<SectionInfo> m_sections;
```

Parse sections in `parsePart`:

```cpp
else if (keyword.find("SOLID SECTION") == 0 || keyword.find("SHELL SECTION") == 0)
{
    std::map<std::string, std::string> sectionParams;
    parseKeywordParams(line, sectionParams);
    
    SectionInfo section;
    section.elementSetName = sectionParams["ELSET"];
    section.materialName = sectionParams["MATERIAL"];
    section.type = (keyword.find("SOLID") == 0) ? "SOLID" : "SHELL";
    
    m_sections.push_back(section);
}
```

Assign materials in `createDomain`:

```cpp
// Find section for this part
for (const auto& section : m_sections)
{
    if (section.elementSetName == someElsetName)
    {
        // Find material by name
        RgMaterial* mat = fem->FindMaterial(section.materialName);
        if (mat)
        {
            domain->SetMaterial(mat);
        }
        break;
    }
}
```

## Testing

### Example Test Case

Create a simple INP file to test:

```
*Heading
Test material conversion
*Material, name=TestMaterial
*Elastic
210000, 0.3
*Plastic
250, 0.0
350, 0.05
*Density
7850
```

Expected output:
```
Converting material: TestMaterial
  Elastic: E=210000, nu=0.3
  Plastic: 2 hardening points
    Initial yield: 250 at eps_p = 0
  Density: 7850
  Created elasto-plastic material: E=210000, nu=0.3, sy=250, H=2000
```

## Error Handling

The converter performs validation:

1. **Missing elastic properties** - Error, cannot create material
2. **Invalid E or nu** - Error with specific bounds checking
3. **Invalid plastic properties** - Falls back to linear elastic
4. **Missing density** - Warning, continues (density is optional)

## Future Enhancements

Potential additions for more complete Abaqus support:

1. **Hyperelastic materials** - Neo-Hookean, Mooney-Rivlin, etc.
2. **Temperature-dependent properties** - Multi-table data
3. **Anisotropic materials** - Engineering constants, orthotropic
4. **Composite materials** - Layered shells
5. **User materials** - UMAT interface (challenging)
6. **Multilinear hardening** - More complex plastic curves
7. **Damage models** - Progressive damage, failure criteria

## Compilation

Add to your CMakeLists.txt or build system:

```cmake
# Add new source files
set(SOURCES
    ${SOURCES}
    AbaqusMaterialConverter.cpp
    # ... other sources
)
```

## Notes

1. **Coordinate systems**: Abaqus and your solver may use different conventions
2. **Units**: Ensure consistency between Abaqus and FEModel units
3. **Material IDs**: Currently based on order in file; may need explicit IDs
4. **Memory management**: Materials are owned by FEModel, not AbaqusImport

## Contact

For questions or issues with material conversion, refer to the FEModel
material documentation and the Abaqus Keywords Reference Manual.
