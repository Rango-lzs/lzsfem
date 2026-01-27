/*********************************************************************
 * \file   AbaqusImport.h
 * \brief  Import Abaqus INP file and convert to FEModel
 *
 * \author
 * \date   January 2026
 *********************************************************************/

#pragma once
#include "datastructure/Vector3d.h"
#include "femcore/fem_export.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

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
    bool load(const char* filename, FEModel* fem);

private:
    struct BoundaryCondition
    {
        std::string nodeSetName;
        int dof;
        double value;
        int step;
    };

    struct ConcentratedLoad
    {
        std::string nodeSetName;
        int dof;
        double magnitude;
        int step;
    };

    struct DistributedLoad
    {
        std::string surfaceName;
        std::string loadType;  // P (pressure), TRVEC (traction)
        double magnitude;
        Vector3d direction;
        int step;
    };

    struct MaterialProperty
    {
        std::string name;
        std::string type;  // ELASTIC, PLASTIC, etc.
        std::map<std::string, std::vector<double>> properties;
    };

    struct StepInfo
    {
        std::string name;
        std::string procedure;  // STATIC, DYNAMIC, etc.
        double timePeriod;
        double initialTimeIncrement;
        double minTimeIncrement;
        double maxTimeIncrement;
    };

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
        int type;           // Element type identifier
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
    bool parseFile(const char* filename);
    bool processData(FEModel* fem);

    // Keyword parsing functions
    bool parseHeading(std::ifstream& file);
    bool parsePart(std::ifstream& file, const std::string& keywordLine);
    bool parseNode(std::ifstream& file, AbaqusPart& part);
    bool parseElement(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine);
    bool parseNset(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine);
    bool parseElset(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine);
    bool parseSurface(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine);
    bool parseMaterial(std::ifstream& file, const std::string& keywordLine);
    bool parseStep(std::ifstream& file, const std::string& keywordLine);
    bool parseBoundary(std::ifstream& file);
    bool parseLoad(std::ifstream& file);
    bool parseCload(std::ifstream& file);
    bool parseDload(std::ifstream& file);

    // Assembly parsing
    bool parseAssembly(std::ifstream& file);
    bool parseInstance(std::ifstream& file, const std::string& keywordLine);
    bool parseEndInstance(std::ifstream& file);
    bool parseEndAssembly(std::ifstream& file);

    // Conversion functions
    RgDomain* createDomain(const AbaqusPart& part, FEMesh* mesh);
    bool createNodes(const AbaqusPart& part, FEMesh* mesh, std::map<int, int>& nodeMap);
    bool createElements(const AbaqusPart& part, RgDomain* domain, const std::map<int, int>& nodeMap);
    bool createNodeSets(const AbaqusPart& part, FEMesh* mesh, const std::map<int, int>& nodeMap);
    bool createElementSets(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain);
    bool createSurfaces(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain);
    bool createBoundaryConditions(FEModel* fem);
    bool createLoads(FEModel* fem);

    // Element type conversion
    int convertElementType(const std::string& abqType);
    int getElementNodeCount(int elemType);

    // Utility functions
    std::string readKeywordLine(std::ifstream& file);
    void parseKeywordParams(const std::string& line, std::map<std::string, std::string>& params);
    std::vector<std::string> splitString(const std::string& str, char delimiter);
    std::string trimString(const std::string& str);
    std::string toUpper(const std::string& str);
    bool isKeyword(const std::string& line);
    bool skipToNextKeyword(std::ifstream& file);

    // Data validation
    bool validatePart(const AbaqusPart& part);
    bool validateNodeId(int nodeId, const std::map<int, int>& nodeMap);

private:
    std::vector<AbaqusPart> m_parts;
    std::vector<AbaqusInstance> m_instances;
    std::vector<MaterialProperty> m_materials;
    std::vector<BoundaryCondition> m_boundaryConditions;
    std::vector<ConcentratedLoad> m_concentratedLoads;
    std::vector<DistributedLoad> m_distributedLoads;
    std::vector<StepInfo> m_steps;

    // Current parsing context
    std::string m_currentPart;
    std::string m_currentInstance;
    std::string m_currentMaterial;
    int m_currentStep;
    bool m_inAssembly;
    bool m_inStep;

    // Global node/element mappings
    std::map<int, int> m_globalNodeMap;  // Abaqus ID -> FEModel index
    std::map<int, int> m_globalElemMap;  // Abaqus ID -> FEModel index

    // Statistics
    int m_totalNodes;
    int m_totalElements;

    // Error handling
    std::string m_lastError;
    int m_lineNumber;

    //int m_globalNodeOffset;  // 全局节点偏移
    //int m_globalElemOffset;  // 全局单元偏移

    private:
    // 全局编号管理
    int m_globalNodeOffset;  // 当前全局节点偏移
    int m_globalElemOffset;  // 当前全局单元偏移

    // 辅助函数
    void applyRotation(Vector3d& pos, const double rotation[7]);

    // 修改后的函数签名
    bool createNodes(const AbaqusPart& part, FEMesh* mesh, std::map<int, int>& nodeMap, const AbaqusInstance& instance);

    bool createElements(const AbaqusPart& part, RgDomain* domain, const std::map<int, int>& nodeMap);

    bool createNodeSets(const AbaqusPart& part, FEMesh* mesh, const std::map<int, int>& nodeMap);

    bool createElementSets(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain);
};