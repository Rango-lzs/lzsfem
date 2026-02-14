/*********************************************************************
 * \file   AbaqusMaterialConverter.h
 * \brief  Convert Abaqus material definitions to FEModel materials
 *
 * \author
 * \date   February 2026
 *********************************************************************/

#pragma once

#include "femcore/fem_export.h"
#include "materials/RgMaterial.h"
#include <string>
#include <vector>
#include <map>

class FEModel;
struct MaterialProperty;

/**
 * @brief Helper class to convert Abaqus material properties to FEModel materials
 */
class FEM_EXPORT AbaqusMaterialConverter
{
public:  
    /**
     * @brief Convert Abaqus material to FEModel material
     * @param matProp Material properties from Abaqus INP
     * @param fem Pointer to FEModel
     * @return Pointer to created RgMaterial, nullptr if conversion failed
     */
    static RgMaterial* ConvertMaterial(const MaterialProperty& matProp, FEModel* fem);

private:
    /**
     * @brief Create linear elastic material
     */
    static RgMaterial* CreateLinearElastic(const MaterialProperty& matProp);

    /**
     * @brief Create elasto-plastic material
     */
    static RgMaterial* CreateElastoPlastic(const MaterialProperty& matProp);

    /**
     * @brief Extract density from material properties
     */
    static double ExtractDensity(const MaterialProperty& matProp);

    /**
     * @brief Check if material has plastic properties
     */
    static bool HasPlasticProperties(const MaterialProperty& matProp);

    /**
     * @brief Get elastic properties (E, nu)
     */
    static bool GetElasticProperties(const MaterialProperty& matProp, 
                                     double& E, double& nu);

    /**
     * @brief Get plastic properties (yield stress, hardening)
     */
    static bool GetPlasticProperties(const MaterialProperty& matProp,
                                     double& sy, double& H);
};
