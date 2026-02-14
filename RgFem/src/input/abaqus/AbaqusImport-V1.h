/*********************************************************************
 * \file   AbaqusImport.h
 * \brief  Import Abaqus INP file and convert to FEModel
 *
 * \author
 * \date   January 2026
 *********************************************************************/

#pragma once
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include "femcore/fem_export.h"

class FEModel;
class FEMesh;
class RgDomain;
class FENode;
class FENodeSet;
class RgElementSet;
class FEFacetSet;

/**
 * @brief Abaqus INP file importer
 *
 * This class reads Abaqus INP files and converts them to FEModel format.
 * Each Abaqus Part corresponds to one RgDomain in the model.
 */
class FEM_EXPORT AbaqusImport
{
public:
    AbaqusImport();
    ~AbaqusImport();

    /**
     * @brief Load an Abaqus INP file and create FEModel
     * @param filename Path to the INP file
     * @param fem Pointer to FEModel to populate
     * @return true if successful, false otherwise
     */
    bool Load(const char* filename, FEModel* fem);

private:
    // Helper structures for parsing
    struct AbaqusNode
    {
        int id;
        double x, y, z;
    };

    struct AbaqusElement
    {
        int id;
        int type;  // Element type identifier
        std::vector<int> nodes;
        std::string elset;  // Element set name
    };

    struct AbaqusPart
    {
        std::string name;
        std::vector<AbaqusNode> nodes;
        std::vector<AbaqusElement> elements;
        std::map<std::string, std::vector<int>> nodeSets;
        std::map<std::string, std::vector<int>> elementSets;
        std::map<std::string, std::vector<std::vector<int>>> surfaces;
    };

    struct AbaqusInstance
    {
        std::string name;
        std::string partName;
        double translation[3];
        double rotation[7];  // axis (3) + angle or rotation matrix
    };

private:
    // Core parsing functions
    bool ParseFile(const char* filename);
    bool ProcessData(FEModel* fem);

    // Keyword parsing functions
    bool ParseHeading(std::ifstream& file);
    bool ParsePart(std::ifstream& file);
    bool ParseNode(std::ifstream& file, AbaqusPart& part);
    bool ParseElement(std::ifstream& file, AbaqusPart& part);
    bool ParseNset(std::ifstream& file, AbaqusPart& part);
    bool ParseElset(std::ifstream& file, AbaqusPart& part);
    bool ParseSurface(std::ifstream& file, AbaqusPart& part);
    bool ParseMaterial(std::ifstream& file);
    bool ParseStep(std::ifstream& file);
    bool ParseBoundary(std::ifstream& file);
    bool ParseLoad(std::ifstream& file);

    // Assembly parsing
    bool ParseAssembly(std::ifstream& file);
    bool ParseInstance(std::ifstream& file);

    // Conversion functions
    RgDomain* CreateDomain(const AbaqusPart& part, FEMesh* mesh);
    bool CreateNodes(const AbaqusPart& part, FEMesh* mesh,
        std::map<int, int>& nodeMap);
    bool CreateElements(const AbaqusPart& part, RgDomain* domain,
        const std::map<int, int>& nodeMap);
    bool CreateNodeSets(const AbaqusPart& part, FEMesh* mesh,
        const std::map<int, int>& nodeMap);
    bool CreateElementSets(const AbaqusPart& part, FEMesh* mesh,
        RgDomain* domain);
    bool CreateSurfaces(const AbaqusPart& part, FEMesh* mesh,
        RgDomain* domain);

    // Element type conversion
    int ConvertElementType(const std::string& abqType);
    int GetElementNodeCount(int elemType);

    // Utility functions
    std::string ReadKeywordLine(std::ifstream& file);
    void ParseKeywordParams(const std::string& line,
        std::map<std::string, std::string>& params);
    std::vector<std::string> SplitString(const std::string& str,
        char delimiter);
    std::string TrimString(const std::string& str);
    std::string ToUpper(const std::string& str);
    bool IsKeyword(const std::string& line);
    bool SkipToNextKeyword(std::ifstream& file);

    // Data validation
    bool ValidatePart(const AbaqusPart& part);
    bool ValidateNodeId(int nodeId, const std::map<int, int>& nodeMap);

private:
    std::vector<AbaqusPart> m_parts;
    std::vector<AbaqusInstance> m_instances;
    std::map<std::string, std::string> m_materials;

    // Current parsing context
    std::string m_currentPart;
    std::string m_currentInstance;
    bool m_inAssembly;

    // Global node/element mappings
    std::map<int, int> m_globalNodeMap;  // Abaqus ID -> FEModel index
    std::map<int, int> m_globalElemMap;  // Abaqus ID -> FEModel index

    // Statistics
    int m_totalNodes;
    int m_totalElements;

    // Error handling
    std::string m_lastError;
    int m_lineNumber;
};