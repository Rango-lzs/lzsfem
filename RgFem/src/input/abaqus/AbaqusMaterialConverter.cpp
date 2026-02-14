/*********************************************************************
 * \file   AbaqusMaterialConverter.cpp
 * \brief  Implementation of Abaqus material converter
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#include "AbaqusMaterialConverter.h"
#include "materials/SmallDeformation/RgLinearElasticMaterial.h"
#include "materials/SmallDeformation/RgElastoPlasticMaterial.h"
#include "femcore/FEModel.h"
#include "logger/log.h"

#include <algorithm>
#include "AbaqusImport.h"

//-----------------------------------------------------------------------------
RgMaterial* AbaqusMaterialConverter::ConvertMaterial(
    const MaterialProperty& matProp, 
    FEModel* fem)
{
    if (!fem)
    {
        RgLogError("FEModel pointer is null");
        return nullptr;
    }

    RgLog("Converting material: %s\n", matProp.name.c_str());

    // Check if material has elastic properties
    double E = 0.0, nu = 0.0;
    if (!GetElasticProperties(matProp, E, nu))
    {
        RgLogError("Material '%s' missing elastic properties", matProp.name.c_str());
        return nullptr;
    }

    RgMaterial* material = nullptr;

    // Check if material has plastic properties
    if (HasPlasticProperties(matProp))
    {
        // Create elasto-plastic material
        double sy = 0.0, H = 0.0;
        if (GetPlasticProperties(matProp, sy, H))
        {
            material = CreateElastoPlastic(matProp);
            RgLog("  Created elasto-plastic material: E=%g, nu=%g, sy=%g, H=%g\n", 
                  E, nu, sy, H);
        }
        else
        {
            RgLogWarning("Plastic properties incomplete for '%s', creating linear elastic instead", 
                        matProp.name.c_str());
            material = CreateLinearElastic(matProp);
        }
    }
    else
    {
        // Create linear elastic material
        material = CreateLinearElastic(matProp);
        RgLog("  Created linear elastic material: E=%g, nu=%g\n", E, nu);
    }

    // Set density if available
    if (material)
    {
        double density = ExtractDensity(matProp);
        if (density > 0.0)
        {
            // Assuming RgMaterial has setDensity method
            // If not, you may need to add it or handle differently
            // material->setDensity(density);
            RgLog("  Density: %g\n", density);
        }
    }

    return material;
}

//-----------------------------------------------------------------------------
RgMaterial* AbaqusMaterialConverter::CreateLinearElastic(const MaterialProperty& matProp)
{
    double E = 0.0, nu = 0.0;
    
    if (!GetElasticProperties(matProp, E, nu))
    {
        return nullptr;
    }

    // Validate material parameters
    if (E <= 0.0)
    {
        RgLogError("Invalid Young's modulus: %g", E);
        return nullptr;
    }

    if (nu < -1.0 || nu >= 0.5)
    {
        RgLogError("Invalid Poisson's ratio: %g (must be in range [-1, 0.5))", nu);
        return nullptr;
    }

    return new SmallDef::RgLinearElastic(E, nu);
}

//-----------------------------------------------------------------------------
RgMaterial* AbaqusMaterialConverter::CreateElastoPlastic(const MaterialProperty& matProp)
{
    double E = 0.0, nu = 0.0, sy = 0.0, H = 0.0;
    
    if (!GetElasticProperties(matProp, E, nu))
    {
        return nullptr;
    }

    if (!GetPlasticProperties(matProp, sy, H))
    {
        return nullptr;
    }

    // Validate material parameters
    if (E <= 0.0)
    {
        RgLogError("Invalid Young's modulus: %g", E);
        return nullptr;
    }

    if (nu < -1.0 || nu >= 0.5)
    {
        RgLogError("Invalid Poisson's ratio: %g", nu);
        return nullptr;
    }

    if (sy <= 0.0)
    {
        RgLogError("Invalid yield stress: %g", sy);
        return nullptr;
    }

    if (H < 0.0)
    {
        RgLogError("Invalid hardening modulus: %g", H);
        return nullptr;
    }

    return new SmallDef::RgElastoPlastic(E, nu, sy, H);
}

//-----------------------------------------------------------------------------
double AbaqusMaterialConverter::ExtractDensity(const MaterialProperty& matProp)
{
    auto it = matProp.properties.find("DENSITY");
    if (it != matProp.properties.end() && !it->second.empty())
    {
        return it->second[0];
    }
    return 0.0;
}

//-----------------------------------------------------------------------------
bool AbaqusMaterialConverter::HasPlasticProperties(const MaterialProperty& matProp)
{
    return (matProp.properties.find("PLASTIC") != matProp.properties.end());
}

//-----------------------------------------------------------------------------
bool AbaqusMaterialConverter::GetElasticProperties(
    const MaterialProperty& matProp, 
    double& E, 
    double& nu)
{
    auto it = matProp.properties.find("ELASTIC");
    if (it != matProp.properties.end())
    {
        const std::vector<double>& props = it->second;
        if (props.size() >= 2)
        {
            E = props[0];   // Young's modulus
            nu = props[1];  // Poisson's ratio
            return true;
        }
        else
        {
            RgLogError("Elastic properties incomplete (need E and nu)");
            return false;
        }
    }
    
    RgLogError("No elastic properties found");
    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusMaterialConverter::GetPlasticProperties(
    const MaterialProperty& matProp,
    double& sy, 
    double& H)
{
    auto it = matProp.properties.find("PLASTIC");
    if (it != matProp.properties.end())
    {
        const std::vector<double>& props = it->second;
        
        // Abaqus PLASTIC format: yield_stress, plastic_strain
        // For first point: yield_stress at zero plastic strain
        if (props.size() >= 1)
        {
            sy = props[0];  // Initial yield stress
            
            // Calculate hardening modulus if multiple points are given
            H = 0.0;  // Default: perfectly plastic
            
            if (props.size() >= 4)  // At least 2 data points
            {
                // Second point: stress, plastic strain
                double stress2 = props[2];
                double pstrain2 = props[3];
                
                // Linear hardening: H = (stress2 - sy) / pstrain2
                if (pstrain2 > 1e-12)
                {
                    H = (stress2 - sy) / pstrain2;
                }
            }
            
            return true;
        }
    }
    
    return false;
}
