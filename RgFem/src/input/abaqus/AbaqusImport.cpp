#include "AbaqusImport.h"

#include "AbaqusBoundaryParser.h"
#include "AbaqusLoadParser.h"
#include "elements/RgElement/RgHex8Element.h"
#include "elements/RgElementSet.h"
#include "femcore/BoundaryCondition/RgBoundaryCondition.h"
#include "femcore/Domain/RgDomain.h"
#include "femcore/Domain/RgSolidDomain.h"
#include "femcore/FEFacetSet.h"
#include "femcore/FEMesh.h"
#include "femcore/FEModel.h"
#include "femcore/FEModelLoad.h"
#include "femcore/FENode.h"
#include "femcore/FENodeSet.h"
#include "femcore/FESolidAnalysis.h"
#include "femcore/Load/RgLoad.h"
#include "femcore/RgLoadController.h"
#include "logger/log.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include "AbaqusMaterialConverter.h"

// ============================================================================
// 单元类型常量定义（如果还没有，需要在某个头文件中定义）
// ============================================================================

enum AbaElementType
{
    FE_HEX8 = 0,
    FE_HEX20,
    FE_HEX27,
    FE_TET4,
    FE_TET10,
    FE_PENTA6,
    FE_PENTA15,
    FE_QUAD4,
    FE_QUAD8,
    FE_TRI3,
    FE_TRI6,
    FE_LINE2,
    // ... 其他类型
};

//-----------------------------------------------------------------------------
AbaqusImport::AbaqusImport()
    : m_inAssembly(false)
    , m_inStep(false)
    , m_currentStep(-1)
    , m_totalNodes(0)
    , m_totalElements(0)
    , m_lineNumber(0)
{
}

//-----------------------------------------------------------------------------
AbaqusImport::~AbaqusImport()
{
}

//-----------------------------------------------------------------------------
bool AbaqusImport::load(const char* filename, FEModel* fem)
{
    if (!filename || !fem)
    {
        m_lastError = "Invalid input parameters";
        // RgLogError(m_lastError.c_str());
        return false;
    }

    m_model = fem;

    // RgLog("Reading Abaqus INP file: %s\n", filename);

    // Clear previous data
    m_parts.clear();
    m_instances.clear();
    m_materials.clear();
    m_boundaryConditions.clear();
    m_concentratedLoads.clear();
    m_distributedLoads.clear();
    m_steps.clear();
    m_globalNodeMap.clear();
    m_globalElemMap.clear();
    m_totalNodes = 0;
    m_totalElements = 0;

    // Parse the INP file
    if (!parseFile(filename))
    {
        RgLogError("Failed to parse INP file: %s", m_lastError.c_str());
        return false;
    }

    // RgLog("Successfully imported %d nodes and %d elements from %d parts\n", m_totalNodes, m_totalElements,
    //       (int)m_parts.size());

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseFile(const char* filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        m_lastError = "Cannot open file";
        return false;
    }

    m_lineNumber = 0;
    std::string line;

    while (std::getline(file, line))
    {
        m_lineNumber++;

        // Skip empty lines and comments
        line = trimString(line);
        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
            continue;

        // Check if this is a keyword line
        if (line[0] == '*' && line[1] != '*')
        {
            std::string keyword = toUpper(line.substr(1));

            // Parse based on keyword
            if (keyword.find("HEADING") == 0)
            {
                parseHeading(file);
            }
            else if (keyword.find("PART") == 0)
            {
                if (!parsePart(file, line))
                    return false;
            }
            else if (keyword.find("ASSEMBLY") == 0)
            {
                m_inAssembly = true;
                if (!parseAssembly(file))
                    return false;

                // Convert Instance data to FEModel
                if (!processData(m_model))
                {
                    RgLogError("Failed to process INP data: %s", m_lastError.c_str());
                    return false;
                }

                m_inAssembly = false;
            }
            else if (keyword.find("MATERIAL") == 0)
            {
                parseMaterial(file, line);
            }
            else if (keyword.find("*BOUNDARY") == 0)
            {
                parseBoundary(file);
            }
            else if (keyword.find("*CLOAD") == 0)
            {
                parseCload(file);
            }
            else if (keyword.find("*DLOAD") == 0)
            {
                parseDload(file);
            }
            else if (keyword.find("STEP") == 0)
            {
                parseStep(file, line);
            }
        }
    }

    file.close();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parsePart(std::ifstream& file, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    AbaqusPart part;
    part.name = params["NAME"];

    if (part.name.empty())
    {
        m_lastError = "Part name is required";
        return false;
    }

    m_currentPart = part.name;
    // RgLog("  Parsing part: %s\n", part.name.c_str());

    std::string line;

    // Parse part contents
    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
            continue;

        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            if (keyword.find("NODE") == 0 && keyword.find("NODE OUTPUT") != 0 && keyword.find("NODE PRINT") != 0)
            {
                if (!parseNode(file, part))
                    return false;
            }
            else if (keyword.find("ELEMENT") == 0 && keyword.find("ELEMENT OUTPUT") != 0)
            {
                if (!parseElement(file, part, line))
                    return false;
            }
            else if (keyword.find("NSET") == 0)
            {
                if (!parseNset(file, part, line))
                    return false;
            }
            else if (keyword.find("ELSET") == 0)
            {
                if (!parseElset(file, part, line))
                    return false;
            }
            else if (keyword.find("SURFACE") == 0)
            {
                if (!parseSurface(file, part, line))
                    return false;
            }
            else if (keyword.find("SOLID SECTION") == 0 || keyword.find("SHELL SECTION") == 0)
            {
                // Parse section definition (maps element set to material)
                std::map<std::string, std::string> sectionParams;
                parseKeywordParams(line, sectionParams);

                std::string elsetName = sectionParams["ELSET"];
                std::string materialName = sectionParams["MATERIAL"];

                // RgLog("    Section: ELSET=%s, MATERIAL=%s\n", elsetName.c_str(), materialName.c_str());

                // Store section information (you may want to create a dedicated structure)
                // For now, we'll skip reading section data lines
                std::string dataLine;
                std::streampos lastPos = file.tellg();
                if (std::getline(file, dataLine))
                {
                    m_lineNumber++;
                    dataLine = trimString(dataLine);
                    if (dataLine.empty() || dataLine[0] == '*')
                    {
                        // Next keyword found, put it back
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                    // else: section data line, can be processed if needed
                }
            }
            else if (keyword.find("END PART") == 0)
            {
                break;
            }
        }
    }

    if (validatePart(part))
    {
        m_parts.push_back(part);
        m_currentPart.clear();
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseNode(std::ifstream& file, AbaqusPart& part)
{
    std::string line;
    int nodeCount = 0;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            // Put the line back for next keyword
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse node data: ID, X, Y, Z
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() < 4)
        {
            m_lastError = "Invalid node data format";
            return false;
        }

        AbaqusNode node;
        node.id = std::stoi(trimString(tokens[0]));
        node.x = std::stod(trimString(tokens[1]));
        node.y = std::stod(trimString(tokens[2]));
        node.z = std::stod(trimString(tokens[3]));

        part.nodes.push_back(node);
        nodeCount++;
        lastPos = file.tellg();
    }

    // RgLog("    Read %d nodes\n", nodeCount);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseElement(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string elemType = params["TYPE"];
    std::string elsetName = params["ELSET"];

    if (elemType.empty())
    {
        m_lastError = "Element type is required";
        return false;
    }

    int typeId = convertElementType(elemType);
    int nodeCount = getElementNodeCount(typeId);
    int elemCount = 0;

    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse element connectivity
        std::vector<std::string> tokens = splitString(line, ',');

        AbaqusElement elem;
        elem.id = std::stoi(trimString(tokens[0]));
        elem.type = typeId;
        elem.elset = elsetName;

        // Read node connectivity
        for (size_t i = 1; i < tokens.size(); i++)
        {
            elem.nodes.push_back(std::stoi(trimString(tokens[i])));
        }

        // Check if element spans multiple lines
        lastPos = file.tellg();
        while (elem.nodes.size() < nodeCount && std::getline(file, line))
        {
            m_lineNumber++;
            line = trimString(line);
            if (line.empty())
            {
                lastPos = file.tellg();
                continue;
            }
            if (line[0] == '*')
            {
                file.seekg(lastPos);
                m_lineNumber--;
                break;
            }

            tokens = splitString(line, ',');
            for (const auto& token : tokens)
            {
                std::string trimmed = trimString(token);
                if (!trimmed.empty())
                    elem.nodes.push_back(std::stoi(trimmed));
            }
            lastPos = file.tellg();
        }

        part.elements.push_back(elem);
        elemCount++;
        lastPos = file.tellg();
    }

    // RgLog("    Read %d %s elements\n", elemCount, elemType.c_str());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseNset(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string nsetName = params["NSET"];
    if (nsetName.empty())
    {
        m_lastError = "NSET name is required";
        return false;
    }

    std::vector<int> nodeIds;
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse node IDs
        std::vector<std::string> tokens = splitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = trimString(token);
            if (!trimmed.empty())
            {
                nodeIds.push_back(std::stoi(trimmed));
            }
        }
        lastPos = file.tellg();
    }

    part.nodeSets[nsetName] = nodeIds;
    // RgLog("    Read node set '%s' with %d nodes\n", nsetName.c_str(), (int)nodeIds.size());
    return true;
}

bool AbaqusImport::parseNset(std::ifstream& file, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string nsetName = params["NSET"];
    if (nsetName.empty())
    {
        m_lastError = "NSET name is required";
        return false;
    }

    std::string insName = params["INSTANCE"];
    AbaqusInstance* partIns = nullptr;
    for (auto& ins : m_instances)
    {
        if (ins.name == insName)
        {
            partIns = &ins;
            break;
        }
    }
    if (!partIns)
        return false;

    std::vector<int> nodeIds;
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse node IDs
        std::vector<std::string> tokens = splitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = trimString(token);
            if (!trimmed.empty())
            {
                nodeIds.push_back(std::stoi(trimmed));
            }
        }
        lastPos = file.tellg();
    }

    partIns->nodeSets.emplace(nsetName, nodeIds);
    // RgLog("    Read node set '%s' with %d nodes\n", nsetName.c_str(), (int)nodeIds.size());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseElset(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string elsetName = params["ELSET"];
    if (elsetName.empty())
    {
        m_lastError = "ELSET name is required";
        return false;
    }

    std::vector<int> elemIds;
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse element IDs
        std::vector<std::string> tokens = splitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = trimString(token);
            if (!trimmed.empty())
            {
                elemIds.push_back(std::stoi(trimmed));
            }
        }
        lastPos = file.tellg();
    }

    part.elementSets[elsetName] = elemIds;
    // RgLog("    Read element set '%s' with %d elements\n", elsetName.c_str(), (int)elemIds.size());
    return true;
}

bool AbaqusImport::parseElset(std::ifstream& file, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string elsetName = params["ELSET"];
    if (elsetName.empty())
    {
        m_lastError = "ELSET name is required";
        return false;
    }

    std::string insName = params["instance"];
    AbaqusInstance partIns{};
    for (auto ins : m_instances)
    {
        if (ins.name == insName)
        {
            partIns = ins;
        }
    }

    std::vector<int> elemIds;
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse element IDs
        std::vector<std::string> tokens = splitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = trimString(token);
            if (!trimmed.empty())
            {
                elemIds.push_back(std::stoi(trimmed));
            }
        }
        lastPos = file.tellg();
    }

    partIns.elementSets[elsetName] = elemIds;
    // RgLog("    Read element set '%s' with %d elements\n", elsetName.c_str(), (int)elemIds.size());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseSurface(std::ifstream& file, AbaqusPart& part, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    std::string surfaceName = params["NAME"];
    if (surfaceName.empty())
    {
        m_lastError = "Surface name is required";
        return false;
    }

    std::vector<std::vector<int>> surfaceData;
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse surface definition: element_set, face_id
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() >= 2)
        {
            std::vector<int> faceData;
            // First token might be element ID or element set name
            // For simplicity, we'll store as placeholder
            // TODO: Properly handle element set references
            faceData.push_back(0);                                 // Placeholder for element/elset reference
            faceData.push_back(std::stoi(trimString(tokens[1])));  // Face ID
            surfaceData.push_back(faceData);
        }
        lastPos = file.tellg();
    }

    part.surfaces[surfaceName] = surfaceData;
    // RgLog("    Read surface '%s' with %d faces\n", surfaceName.c_str(), (int)surfaceData.size());
    return true;
}

////-----------------------------------------------------------------------------
//bool AbaqusImport::parseMaterial(std::ifstream& file, const std::string& keywordLine)
//{
//    std::map<std::string, std::string> params;
//    parseKeywordParams(keywordLine, params);
//
//    MaterialProperty material;
//    material.name = params["NAME"];
//
//    if (material.name.empty())
//    {
//        m_lastError = "Material name is required";
//        return false;
//    }
//
//    m_currentMaterial = material.name;
//    // RgLog("  Parsing material: %s\n", material.name.c_str());
//
//    std::string line;
//    std::streampos lastPos = file.tellg();
//
//    while (std::getline(file, line))
//    {
//        m_lineNumber++;
//        line = trimString(line);
//
//        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
//        {
//            lastPos = file.tellg();
//            continue;
//        }
//
//        if (line[0] == '*')
//        {
//            std::string keyword = toUpper(line.substr(1));
//
//            if (keyword.find("ELASTIC") == 0)
//            {
//                // Read elastic properties
//                lastPos = file.tellg();
//                if (std::getline(file, line))
//                {
//                    m_lineNumber++;
//                    line = trimString(line);
//                    if (!line.empty() && line[0] != '*')
//                    {
//                        std::vector<std::string> tokens = splitString(line, ',');
//                        if (tokens.size() >= 2)
//                        {
//                            std::vector<double> elasticProps;
//                            elasticProps.push_back(std::stod(trimString(tokens[0])));  // E
//                            elasticProps.push_back(std::stod(trimString(tokens[1])));  // nu
//                            material.properties["ELASTIC"] = elasticProps;
//                            material.type = "ELASTIC";
//                        }
//                        lastPos = file.tellg();
//                    }
//                    else
//                    {
//                        file.seekg(lastPos);
//                        m_lineNumber--;
//                    }
//                }
//            }
//            else if (keyword.find("DENSITY") == 0)
//            {
//                // Read density
//                lastPos = file.tellg();
//                if (std::getline(file, line))
//                {
//                    m_lineNumber++;
//                    line = trimString(line);
//                    if (!line.empty() && line[0] != '*')
//                    {
//                        std::vector<double> densityProp;
//                        densityProp.push_back(std::stod(trimString(line)));
//                        material.properties["DENSITY"] = densityProp;
//                        lastPos = file.tellg();
//                    }
//                    else
//                    {
//                        file.seekg(lastPos);
//                        m_lineNumber--;
//                    }
//                }
//            }
//            else
//            {
//                // Unknown material property, return to process next material/keyword
//                file.seekg(lastPos);
//                m_lineNumber--;
//                break;
//            }
//        }
//    }
//
//    m_materials.push_back(material);
//    m_currentMaterial.clear();
//    return true;
//}

bool AbaqusImport::parseMaterial(std::ifstream& file, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    MaterialProperty material;
    material.name = params["NAME"];

    if (material.name.empty())
    {
        m_lastError = "Material name is required";
        return false;
    }

    m_currentMaterial = material.name;
    RgLog("  Parsing material: %s\n", material.name.c_str());

    std::string line;
    std::streampos lastPos = file.tellg();

    // Main parsing loop for this material
    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        // Skip empty lines and comments
        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
        {
            lastPos = file.tellg();
            continue;
        }

        // Check for keyword line
        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            // === ELASTIC PROPERTIES ===
            if (keyword.find("ELASTIC") == 0)
            {
                // Read elastic properties: E, nu
                lastPos = file.tellg();
                if (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);
                    if (!line.empty() && line[0] != '*')
                    {
                        std::vector<std::string> tokens = splitString(line, ',');
                        if (tokens.size() >= 2)
                        {
                            std::vector<double> elasticProps;
                            elasticProps.push_back(std::stod(trimString(tokens[0])));  // E
                            elasticProps.push_back(std::stod(trimString(tokens[1])));  // nu
                            material.properties["ELASTIC"] = elasticProps;

                            // Only set type if not already set (plastic takes precedence)
                            if (material.type.empty())
                            {
                                material.type = "ELASTIC";
                            }

                            RgLog("    Elastic: E=%s, nu=%s\n", trimString(tokens[0]).c_str(),
                                  trimString(tokens[1]).c_str());
                        }
                        lastPos = file.tellg();
                    }
                    else
                    {
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                }
            }
            // === PLASTIC PROPERTIES ===
            else if (keyword.find("PLASTIC") == 0)
            {
                // Read plastic properties (hardening curve)
                // Format: yield_stress, plastic_strain (one or more lines)
                std::vector<double> plasticProps;

                lastPos = file.tellg();
                while (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);

                    if (line.empty())
                    {
                        lastPos = file.tellg();
                        continue;
                    }

                    if (line[0] == '*')
                    {
                        // Next keyword, put it back and break
                        file.seekg(lastPos);
                        m_lineNumber--;
                        break;
                    }

                    // Parse plastic data: yield_stress, plastic_strain
                    std::vector<std::string> tokens = splitString(line, ',');
                    if (tokens.size() >= 2)
                    {
                        plasticProps.push_back(std::stod(trimString(tokens[0])));  // stress
                        plasticProps.push_back(std::stod(trimString(tokens[1])));  // plastic strain
                    }
                    else if (tokens.size() == 1)
                    {
                        // Only yield stress given (perfectly plastic)
                        plasticProps.push_back(std::stod(trimString(tokens[0])));
                        plasticProps.push_back(0.0);  // zero plastic strain
                    }

                    lastPos = file.tellg();
                }

                if (!plasticProps.empty())
                {
                    material.properties["PLASTIC"] = plasticProps;
                    material.type = "PLASTIC";

                    RgLog("    Plastic: %d hardening points\n", (int)(plasticProps.size() / 2));
                    if (plasticProps.size() >= 2)
                    {
                        RgLog("      Initial yield: %g at eps_p = %g\n", plasticProps[0], plasticProps[1]);
                    }
                }
            }
            // === DENSITY ===
            else if (keyword.find("DENSITY") == 0)
            {
                // Read density
                lastPos = file.tellg();
                if (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);
                    if (!line.empty() && line[0] != '*')
                    {
                        std::vector<double> densityProp;
                        double rho = std::stod(trimString(line));
                        densityProp.push_back(rho);
                        material.properties["DENSITY"] = densityProp;

                        RgLog("    Density: %g\n", rho);
                        lastPos = file.tellg();
                    }
                    else
                    {
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                }
            }
            // === CHECK FOR EXIT CONDITIONS ===
            else
            {
                // Check if this is a major keyword that indicates end of material definition
                bool shouldExit = false;

                // List of keywords that end the material definition
                if (keyword.find("MATERIAL") == 0)  // Next material
                {
                    shouldExit = true;
                }
                else if (keyword.find("PART") == 0)  // Part definition
                {
                    shouldExit = true;
                }
                else if (keyword.find("ASSEMBLY") == 0)  // Assembly section
                {
                    shouldExit = true;
                }
                else if (keyword.find("STEP") == 0)  // Analysis step
                {
                    shouldExit = true;
                }
                else if (keyword.find("BOUNDARY") == 0)  // Boundary conditions
                {
                    shouldExit = true;
                }
                else if (keyword.find("CLOAD") == 0 || keyword.find("DLOAD") == 0)  // Loads
                {
                    shouldExit = true;
                }
                else if (keyword.find("OUTPUT") == 0 || keyword.find("NODE OUTPUT") == 0 ||
                         keyword.find("ELEMENT OUTPUT") == 0)  // Output requests
                {
                    shouldExit = true;
                }
                else if (keyword.find("SECTION") == 0)  // Section definitions
                {
                    shouldExit = true;
                }
                else if (keyword.find("HEADING") == 0)  // Heading
                {
                    shouldExit = true;
                }
                else if (keyword.find("PREPRINT") == 0)  // Preprint
                {
                    shouldExit = true;
                }

                if (shouldExit)
                {
                    // Major keyword found - put it back and exit
                    file.seekg(lastPos);
                    m_lineNumber--;
                    RgLog("  Material '%s' definition ended at keyword: *%s\n", material.name.c_str(), keyword.c_str());
                    break;
                }

                // Unknown/unsupported material property keyword
                // Skip it and its data lines
                RgLogWarning("Unsupported material property: *%s (line %d) - skipping", keyword.c_str(), m_lineNumber);

                // Skip all data lines for this unknown keyword
                lastPos = file.tellg();
                while (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);

                    if (line.empty())
                    {
                        lastPos = file.tellg();
                        continue;
                    }

                    if (line[0] == '*')
                    {
                        // Hit next keyword, put it back and break
                        file.seekg(lastPos);
                        m_lineNumber--;
                        break;
                    }

                    // Data line for unknown keyword - just skip it
                    RgLog("      Skipping data: %s\n", line.c_str());
                    lastPos = file.tellg();
                }
            }
        }
    }

    // Validate material has at least elastic properties
    if (material.properties.find("ELASTIC") == material.properties.end())
    {
        RgLogWarning("Material '%s' has no elastic properties defined", material.name.c_str());
    }

    m_materials.push_back(material);
    m_currentMaterial.clear();

    RgLog("  Material '%s' parsed successfully (%s)\n", material.name.c_str(),
          material.type.empty() ? "UNDEFINED" : material.type.c_str());

    return true;
}

// ============================================================================
// NEW function to add to AbaqusImport: Create FEModel materials
// ============================================================================

bool AbaqusImport::createMaterials(FEModel* fem, const MaterialProperty& matProp)
{
    RgLog("\n");
    RgLog("=================================================\n");
    RgLog("Creating Materials\n");
    RgLog("=================================================\n");
    RgLog("Converting material '%s'...\n", matProp.name.c_str());

    // Use the converter to create FEModel material
    RgMaterial* material = AbaqusMaterialConverter::ConvertMaterial(matProp, fem);
    if (material)
    {
        // Set material ID based on order (1-based indexing)
        // Note: Assuming RgMaterial or FEMaterial has SetID method
        // If not available, you may need to add it or handle differently
        // material->SetID(i + 1);

        // Note: Assuming RgMaterial has SetName method
        // You might need to add this to your RgMaterial base class
        // material->SetName(matProp.name);

        // Add material to model
        fem->AddMaterial(material);
        RgLog("    Successfully created\n");
    }
    else
    {
        RgLogError("Failed to convert material '%s'", matProp.name.c_str());
    }

    return true;
}


//-----------------------------------------------------------------------------
bool AbaqusImport::parseStep(std::ifstream& file, const std::string& keywordLine)
{
    RgLog("=================================================\n");
    RgLog("Parsing Abaqus Step Section\n");
    RgLog("=================================================\n");

    // Parse step keyword parameters
    std::map<std::string, std::string> stepParams;
    parseKeywordParams(keywordLine, stepParams);

    // Create step info structure for storing parsed data
    StepInfo stepInfo;
    stepInfo.name = stepParams["NAME"];
    if (stepInfo.name.empty())
    {
        stepInfo.name = "Step-" + std::to_string(m_currentStep + 1);
    }

    // Default values
    stepInfo.procedure = "STATIC";
    stepInfo.timePeriod = 1.0;
    stepInfo.initialTimeIncrement = 0.1;
    stepInfo.minTimeIncrement = 1e-5;
    stepInfo.maxTimeIncrement = 0.1;

    RgLog("  Step name: %s\n", stepInfo.name.c_str());
    RgLog("  Step number: %d\n", m_currentStep);

    // Create the FEAnalysis object
    // Note: You may need to create specific analysis types based on your implementation
    // For now, we use the base FEAnalysis class
    FEAnalysis* analysisStep = new FESolidAnalysis();
    analysisStep->SetName(stepInfo.name);

    // Temporary storage for BCs and loads defined in this step
    std::vector<RgBoundaryCondition*> stepBCs;
    std::vector<RgLoad*> stepLoads;

    std::string line;
    std::streampos lastPos = file.tellg();

    m_inStep = true;
    bool stepConfigured = false;

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        // Skip empty lines and comments
        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
        {
            lastPos = file.tellg();
            continue;
        }

        // Check for keyword
        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            // Parse step procedure type
            if (keyword.find("STATIC") == 0)
            {
                stepInfo.procedure = "STATIC";
                analysisStep->m_nanalysis = 0;  // STATIC
                RgLog("  Procedure: STATIC\n");

                // Parse STATIC parameters
                std::map<std::string, std::string> procParams;
                parseKeywordParams(line, procParams);

                // Read step control data (next line after *Static)
                lastPos = file.tellg();
                if (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);

                    if (!line.empty() && line[0] != '*')
                    {
                        std::vector<std::string> tokens = splitString(line, ',');

                        // Format: initial_increment, period, min_increment, max_increment
                        if (tokens.size() >= 1)
                        {
                            stepInfo.initialTimeIncrement = std::stod(trimString(tokens[0]));
                            analysisStep->m_dt0 = stepInfo.initialTimeIncrement;
                            RgLog("    Initial time increment: %g\n", stepInfo.initialTimeIncrement);
                        }
                        if (tokens.size() >= 2)
                        {
                            stepInfo.timePeriod = std::stod(trimString(tokens[1]));
                            analysisStep->m_final_time = stepInfo.timePeriod;
                            RgLog("    Time period: %g\n", stepInfo.timePeriod);
                        }
                        if (tokens.size() >= 3)
                        {
                            stepInfo.minTimeIncrement = std::stod(trimString(tokens[2]));
                            RgLog("    Min time increment: %g\n", stepInfo.minTimeIncrement);
                        }
                        if (tokens.size() >= 4)
                        {
                            stepInfo.maxTimeIncrement = std::stod(trimString(tokens[3]));
                            RgLog("    Max time increment: %g\n", stepInfo.maxTimeIncrement);
                        }

                        stepConfigured = true;
                        lastPos = file.tellg();
                    }
                    else
                    {
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                }
            }
            else if (keyword.find("DYNAMIC") == 0)
            {
                stepInfo.procedure = "DYNAMIC";
                analysisStep->m_nanalysis = 1;  // DYNAMIC
                RgLog("  Procedure: DYNAMIC\n");

                // Parse DYNAMIC parameters
                std::map<std::string, std::string> procParams;
                parseKeywordParams(line, procParams);

                // Read dynamic step parameters
                lastPos = file.tellg();
                if (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);

                    if (!line.empty() && line[0] != '*')
                    {
                        std::vector<std::string> tokens = splitString(line, ',');

                        if (tokens.size() >= 1)
                        {
                            stepInfo.initialTimeIncrement = std::stod(trimString(tokens[0]));
                            analysisStep->m_dt0 = stepInfo.initialTimeIncrement;
                            RgLog("    Initial time increment: %g\n", stepInfo.initialTimeIncrement);
                        }
                        if (tokens.size() >= 2)
                        {
                            stepInfo.timePeriod = std::stod(trimString(tokens[1]));
                            analysisStep->m_final_time = stepInfo.timePeriod;
                            RgLog("    Time period: %g\n", stepInfo.timePeriod);
                        }

                        stepConfigured = true;
                        lastPos = file.tellg();
                    }
                    else
                    {
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                }
            }
            else if (keyword.find("BOUNDARY") == 0)
            {
                RgLog("  Parsing BOUNDARY conditions in step\n");

                // Store current BC count in model
                int bcCountBefore = m_model->BoundaryConditions();

                // Parse boundary conditions
                // Note: parseBoundary() should add BCs directly to m_model
                if (!parseBoundary(file))
                {
                    RgLogWarning("    Failed to parse boundary conditions\n");
                }

                // Collect newly added BCs
                int bcCountAfter = m_model->BoundaryConditions();
                for (int i = bcCountBefore; i < bcCountAfter; ++i)
                {
                    RgBoundaryCondition* bc = m_model->BoundaryCondition(i);
                    stepBCs.push_back(bc);
                    RgLog("    Added BC '%s' to step\n", bc->GetName().c_str());
                }
            }
            else if (keyword.find("CLOAD") == 0)
            {
                RgLog("  Parsing CLOAD in step\n");

                // Store current load count in model
                int loadCountBefore = m_model->ModelLoads();

                // Parse concentrated loads
                if (!parseCload(file))
                {
                    RgLogWarning("    Failed to parse concentrated loads\n");
                }

                // Collect newly added loads
                int loadCountAfter = m_model->ModelLoads();
                for (int i = loadCountBefore; i < loadCountAfter; ++i)
                {
                    RgLoad* load = m_model->ModelLoad(i);
                    stepLoads.push_back(load);
                    RgLog("    Added load '%s' to step\n", load->GetName().c_str());
                }
            }
            else if (keyword.find("DLOAD") == 0)
            {
                RgLog("  Parsing DLOAD in step\n");

                // Store current load count in model
                int loadCountBefore = m_model->ModelLoads();

                // Parse distributed loads
                if (!parseDload(file))
                {
                    RgLogWarning("    Failed to parse distributed loads\n");
                }

                // Collect newly added loads
                int loadCountAfter = m_model->ModelLoads();
                for (int i = loadCountBefore; i < loadCountAfter; ++i)
                {
                    RgLoad* load = m_model->ModelLoad(i);
                    stepLoads.push_back(load);
                    RgLog("    Added load '%s' to step\n", load->GetName().c_str());
                }
            }
            else if (keyword.find("OUTPUT") == 0 || keyword.find("NODE OUTPUT") == 0 ||
                     keyword.find("ELEMENT OUTPUT") == 0 || keyword.find("CONTACT OUTPUT") == 0)
            {
                RgLog("  Found OUTPUT request (skipping)\n");
                skipToNextKeyword(file);
            }
            else if (keyword.find("END STEP") == 0)
            {
                RgLog("  End of step section\n");
                break;
            }
            else
            {
                RgLogWarning("  Unknown keyword in step: %s (skipping)\n", keyword.c_str());
                skipToNextKeyword(file);
            }

            lastPos = file.tellg();
        }
    }

    // Set default time parameters if not configured
    if (!stepConfigured)
    {
        analysisStep->m_dt0 = stepInfo.initialTimeIncrement;
        analysisStep->m_final_time = stepInfo.timePeriod;
        RgLog("  Using default time parameters\n");
    }

    // Calculate number of time steps
    if (analysisStep->m_dt0 > 0)
    {
        analysisStep->m_ntime = (int)(analysisStep->m_final_time / analysisStep->m_dt0);
    }
    else
    {
        analysisStep->m_ntime = 10;  // default
    }

    // Add BCs to the step
    for (auto bc : stepBCs)
    {
        analysisStep->AddBoundaryCondition(bc);
    }

    // Add loads to the step
    for (auto load : stepLoads)
    {
        analysisStep->AddLoad(load);
    }

    // Set up inheritance from previous step
    if (m_currentStep > 0 && m_model->Steps() > 0)
    {
        // Get previous step
        FEAnalysis* prevStep = m_model->GetStep(m_currentStep - 1);
        if (prevStep)
        {
            analysisStep->SetPreviousStep(prevStep);
            analysisStep->SetActivationMode(StepActivationMode::INHERITED);
            RgLog("  Will inherit BCs and loads from previous step\n");
        }
    }
    else
    {
        // First step - no inheritance
        analysisStep->SetActivationMode(StepActivationMode::NEW);
        RgLog("  First step - no inheritance\n");
    }

    // Add step to model
    m_model->AddStep(analysisStep);

    // Store step info for reference
    m_steps.push_back(stepInfo);
    m_currentStep++;
    m_inStep = false;

    // Summary
    RgLog("-------------------------------------------------\n");
    RgLog("Step Summary:\n");
    RgLog("  Name: %s\n", analysisStep->GetName().c_str());
    RgLog("  Type: %s\n", stepInfo.procedure.c_str());
    RgLog("  BCs defined in this step: %d\n", (int)stepBCs.size());
    RgLog("  Loads defined in this step: %d\n", (int)stepLoads.size());
    RgLog("  Time period: %g\n", stepInfo.timePeriod);
    RgLog("  Initial dt: %g\n", stepInfo.initialTimeIncrement);
    RgLog("  Number of time steps: %d\n", analysisStep->m_ntime);
    RgLog("=================================================\n\n");

    return true;
}

//-----------------------------------------------------------------------------
// Helper function to create FEAnalysis from StepInfo
bool AbaqusImport::createFEAnalysisStep(const StepInfo& stepInfo)
{
    if (!m_model)
    {
        m_lastError = "Model pointer is null";
        return false;
    }

    // Create FEAnalysis based on procedure type
    FEAnalysis* analysis = nullptr;

    if (stepInfo.procedure == "STATIC")
    {
        // Create static analysis step
        // You may need to adjust this based on your actual FEAnalysis class hierarchy
        // For example: analysis = new FEStaticAnalysis(m_model);

        // If FEAnalysis is the base class:
        analysis = new FESolidAnalysis();
        analysis->SetName(stepInfo.name);

        // Set time parameters
        // Note: You'll need to check your FEAnalysis interface for exact method names
        // These are examples based on common FE code patterns:

        // analysis->SetTimeSteps(stepInfo.initialTimeIncrement);
        // analysis->SetFinalTime(stepInfo.timePeriod);
        // analysis->SetMinTimeStep(stepInfo.minTimeIncrement);
        // analysis->SetMaxTimeStep(stepInfo.maxTimeIncrement);

        RgLog("  Created STATIC analysis step: %s\n", stepInfo.name.c_str());
        RgLog("    Time period: %g\n", stepInfo.timePeriod);
        RgLog("    Initial dt: %g\n", stepInfo.initialTimeIncrement);
    }
    else if (stepInfo.procedure == "DYNAMIC")
    {
        // Create dynamic analysis step
        // analysis = new FEDynamicAnalysis(m_model);

        analysis = new FESolidAnalysis();
        analysis->SetName(stepInfo.name);

        RgLog("  Created DYNAMIC analysis step: %s\n", stepInfo.name.c_str());
        RgLog("    Time period: %g\n", stepInfo.timePeriod);
        RgLog("    Initial dt: %g\n", stepInfo.initialTimeIncrement);
    }
    else
    {
        RgLogError("Unsupported analysis procedure: %s\n", stepInfo.procedure.c_str());
        return false;
    }

    if (analysis)
    {
        // Add step to model
        m_model->AddStep(analysis);
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseBoundary(std::ifstream& file)
{
    RgLog("Parsing *Boundary section...\n");

    // Read lines until we hit another keyword or EOF
    std::vector<std::string> lines;

    // The *Boundary keyword line should already be consumed by the caller
    // So we add it manually
    lines.push_back("*Boundary");

    std::string line;
    std::streampos lastValidPos = file.tellg();

    int dataLineCount = 0;

    while (std::getline(file, line))
    {
        // Check if this line is a keyword (starts with * but not **)
        std::string trimmed = line;

        // Trim leading whitespace
        size_t firstNonSpace = trimmed.find_first_not_of(" \t\r\n");
        if (firstNonSpace != std::string::npos)
        {
            trimmed = trimmed.substr(firstNonSpace);
        }
        else
        {
            trimmed.clear();
        }

        // Check for keyword
        bool isKeyword = false;
        if (!trimmed.empty() && trimmed[0] == '*')
        {
            if (trimmed.size() == 1)
            {
                // Just '*' - treat as keyword
                isKeyword = true;
            }
            else if (trimmed[1] != '*')
            {
                // Starts with * but not ** - it's a keyword
                isKeyword = true;
            }
            // else: ** is a comment, not a keyword
        }

        if (isKeyword)
        {
            // We've hit the next keyword, stop reading
            // Seek back to before this line
            file.seekg(lastValidPos);
            break;
        }

        // Count non-comment, non-empty lines
        if (!trimmed.empty() && (trimmed.size() < 2 || trimmed.substr(0, 2) != "**"))
        {
            // Add this line to our collection
            lines.push_back(line);
            dataLineCount++;
        }

        // Save position for potential seek back
        lastValidPos = file.tellg();
    }

    if (dataLineCount == 0)
    {
        RgLogWarning("*Boundary section has no data lines");
        return false;
    }

    // Create parser
    AbaqusBoundaryParser parser(m_model);

    // Parse the boundary conditions
    int linesConsumed = parser.ParseBoundary(lines, 0);

    if (linesConsumed == 0)
    {
        RgLogError("Failed to parse *Boundary section");
        return false;
    }

    // Create FE boundary condition objects
    if (!parser.CreateBoundaryConditions())
    {
        RgLogError("Failed to create boundary condition objects");
        return false;
    }

    // Report results
    const auto& bcData = parser.GetBoundaryData();
    RgLog("  Parsed %d boundary condition entries from %d data lines\n", (int)bcData.size(), dataLineCount);

    // Detailed logging (optional - can be controlled by verbosity level)
#ifdef VERBOSE_LOGGING
    for (size_t i = 0; i < bcData.size(); ++i)
    {
        const auto& data = bcData[i];
        if (data.isEncastre)
        {
            RgLog("    [%d] ENCASTRE: %s\n", (int)i + 1, data.nodeSetName.c_str());
        }
        else
        {
            RgLog("    [%d] %s: DOF %d-%d = %g\n", (int)i + 1, data.nodeSetName.c_str(), data.firstDOF, data.lastDOF,
                  data.magnitude);
        }
    }
#endif

    RgLog("  *Boundary section parsed successfully\n");
}


// Parse *Cload (Concentrated Load)
bool AbaqusImport::parseCload(std::ifstream& file)
{
    RgLog("Parsing *Cload section...\n");

    // Read lines until we hit another keyword or EOF
    std::vector<std::string> lines;

    // Add the *Cload keyword line
    lines.push_back("*Cload");

    std::string line;
    std::streampos lastValidPos = file.tellg();

    int dataLineCount = 0;

    while (std::getline(file, line))
    {
        // Trim and check for keyword
        std::string trimmed = line;
        size_t firstNonSpace = trimmed.find_first_not_of(" \t\r\n");
        if (firstNonSpace != std::string::npos)
        {
            trimmed = trimmed.substr(firstNonSpace);
        }
        else
        {
            trimmed.clear();
        }

        // Check for keyword
        bool isKeyword = false;
        if (!trimmed.empty() && trimmed[0] == '*')
        {
            if (trimmed.size() == 1 || trimmed[1] != '*')
            {
                isKeyword = true;
            }
        }

        if (isKeyword)
        {
            file.seekg(lastValidPos);
            break;
        }

        lines.push_back(line);

        if (!trimmed.empty() && (trimmed.size() < 2 || trimmed.substr(0, 2) != "**"))
        {
            dataLineCount++;
        }

        lastValidPos = file.tellg();
    }

    if (dataLineCount == 0)
    {
        RgLogWarning("*Cload section has no data lines");
        return false;
    }

    // Create parser
    AbaqusLoadParser parser(m_model);

    // Parse concentrated loads
    int linesConsumed = parser.ParseCload(lines, 0);

    if (linesConsumed == 0)
    {
        RgLogError("Failed to parse *Cload section");
        return false;
    }

    // Create load objects
    if (!parser.CreateLoads())
    {
        RgLogError("Failed to create concentrated loads");
        return false;
    }

    // Report results
    const auto& loadData = parser.GetLoadData();
    RgLog("  Parsed %d concentrated load entries from %d data lines\n", (int)loadData.size(), dataLineCount);

#ifdef VERBOSE_LOGGING
    for (size_t i = 0; i < loadData.size(); ++i)
    {
        const auto& data = loadData[i];
        RgLog("    [%d] %s: DOF %d = %g\n", (int)i + 1, data.targetName.c_str(), data.dof, data.magnitude);
    }
#endif

    RgLog("  *Cload section parsed successfully\n");
}

//-----------------------------------------------------------------------------
// Parse *Dload (Distributed Load)
bool AbaqusImport::parseDload(std::ifstream& file)
{
    RgLog("Parsing *Dload section...\n");

    // Read lines until we hit another keyword or EOF
    std::vector<std::string> lines;

    // Add the *Dload keyword line
    lines.push_back("*Dload");

    std::string line;
    std::streampos lastValidPos = file.tellg();

    int dataLineCount = 0;

    while (std::getline(file, line))
    {
        // Trim and check for keyword
        std::string trimmed = line;
        size_t firstNonSpace = trimmed.find_first_not_of(" \t\r\n");
        if (firstNonSpace != std::string::npos)
        {
            trimmed = trimmed.substr(firstNonSpace);
        }
        else
        {
            trimmed.clear();
        }

        // Check for keyword
        bool isKeyword = false;
        if (!trimmed.empty() && trimmed[0] == '*')
        {
            if (trimmed.size() == 1 || trimmed[1] != '*')
            {
                isKeyword = true;
            }
        }

        if (isKeyword)
        {
            file.seekg(lastValidPos);
            break;
        }

        lines.push_back(line);

        if (!trimmed.empty() && (trimmed.size() < 2 || trimmed.substr(0, 2) != "**"))
        {
            dataLineCount++;
        }

        lastValidPos = file.tellg();
    }

    if (dataLineCount == 0)
    {
        RgLogWarning("*Dload section has no data lines");
        return false;
    }

    // Create parser
    AbaqusLoadParser parser(m_model);

    // Parse distributed loads
    int linesConsumed = parser.ParseDload(lines, 0);

    if (linesConsumed == 0)
    {
        RgLogError("Failed to parse *Dload section");
        return false;
    }

    // Create load objects
    if (!parser.CreateLoads())
    {
        RgLogError("Failed to create distributed loads");
        return false;
    }

    // Report results
    const auto& loadData = parser.GetLoadData();
    RgLog("  Parsed %d distributed load entries from %d data lines\n", (int)loadData.size(), dataLineCount);

#ifdef VERBOSE_LOGGING
    for (size_t i = 0; i < loadData.size(); ++i)
    {
        const auto& data = loadData[i];
        RgLog("    [%d] %s on %s: type=%s, magnitude=%g\n", (int)i + 1, data.name.c_str(), data.targetName.c_str(),
              data.loadType.c_str(), data.magnitude);
    }
#endif

    RgLog("  *Dload section parsed successfully\n");
}

//-----------------------------------------------------------------------------
// Generic parseLoad that can handle both *Cload and *Dload
// Call this if you want a unified interface
// bool AbaqusImport::parseLoad(std::ifstream& file, const std::string& loadType)
//{
//    // Convert to uppercase for comparison
//    std::string upperType = loadType;
//    std::transform(upperType.begin(), upperType.end(), upperType.begin(), ::toupper);
//
//    if (upperType.find("CLOAD") != std::string::npos)
//    {
//        parseCload(file);
//    }
//    else if (upperType.find("DLOAD") != std::string::npos)
//    {
//        parseDload(file);
//    }
//    else
//    {
//        RgLogError("Unknown load type: %s", loadType.c_str());
//    }
//}

////-----------------------------------------------------------------------------
// bool AbaqusImport::parseCload(std::ifstream& file)
//{
//     std::string line;
//     std::streampos lastPos = file.tellg();
//
//     while (std::getline(file, line))
//     {
//         m_lineNumber++;
//         line = trimString(line);
//
//         if (line.empty())
//         {
//             lastPos = file.tellg();
//             continue;
//         }
//         if (line[0] == '*')
//         {
//             file.seekg(lastPos);
//             m_lineNumber--;
//             break;
//         }
//
//         // Parse concentrated load: node_set, dof, magnitude
//         std::vector<std::string> tokens = splitString(line, ',');
//         if (tokens.size() >= 3)
//         {
//             ConcentratedLoad load;
//             load.nodeSetName = trimString(tokens[0]);
//             load.dof = std::stoi(trimString(tokens[1]));
//             load.magnitude = std::stod(trimString(tokens[2]));
//             load.step = m_currentStep;
//             m_concentratedLoads.push_back(load);
//         }
//         lastPos = file.tellg();
//     }
//
//     return true;
// }
//
////-----------------------------------------------------------------------------
// bool AbaqusImport::parseDload(std::ifstream& file)
//{
//     std::string line;
//     std::streampos lastPos = file.tellg();
//
//     while (std::getline(file, line))
//     {
//         m_lineNumber++;
//         line = trimString(line);
//
//         if (line.empty())
//         {
//             lastPos = file.tellg();
//             continue;
//         }
//         if (line[0] == '*')
//         {
//             file.seekg(lastPos);
//             m_lineNumber--;
//             break;
//         }
//
//         // Parse distributed load: surface, load_type, magnitude
//         std::vector<std::string> tokens = splitString(line, ',');
//         if (tokens.size() >= 3)
//         {
//             DistributedLoad load;
//             load.surfaceName = trimString(tokens[0]);
//             load.loadType = trimString(tokens[1]);
//             load.magnitude = std::stod(trimString(tokens[2]));
//             load.step = m_currentStep;
//             m_distributedLoads.push_back(load);
//         }
//         lastPos = file.tellg();
//     }
//
//     return true;
// }

//-----------------------------------------------------------------------------
bool AbaqusImport::parseHeading(std::ifstream& file)
{
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // RgLog("Heading: %s\n", line.c_str());
        lastPos = file.tellg();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseAssembly(std::ifstream& file)
{
    std::string line;

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
            continue;

        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            if (keyword.find("INSTANCE") == 0)
            {
                if (!parseInstance(file, line))
                    return false;
            }
            else if (keyword.find("NSET") == 0)
            {
                if (!parseNset(file, line))
                    return false;
            }
            else if (keyword.find("ELSET") == 0)
            {
                if (!parseElset(file, line))
                    return false;
            }
            else if (keyword.find("END ASSEMBLY") == 0)
            {
                break;
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseInstance(std::ifstream& file, const std::string& keywordLine)
{
    std::map<std::string, std::string> params;
    parseKeywordParams(keywordLine, params);

    AbaqusInstance instance;
    instance.name = params["NAME"];
    instance.partName = params["PART"];

    if (instance.name.empty() || instance.partName.empty())
    {
        m_lastError = "Instance must have name and part";
        return false;
    }

    // RgLog("  Parsing instance: %s (part: %s)\n", instance.name.c_str(), instance.partName.c_str());

    // Initialize transformation to identity
    instance.translation[0] = instance.translation[1] = instance.translation[2] = 0.0;
    for (int i = 0; i < 7; i++)
        instance.rotation[i] = 0.0;

    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty())
        {
            lastPos = file.tellg();
            continue;
        }
        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));
            if (keyword.find("END INSTANCE") == 0)
            {
                break;
            }
            file.seekg(lastPos);
            m_lineNumber--;
            break;
        }

        // Parse transformation data (translation, rotation)
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() >= 3)
        {
            instance.translation[0] = std::stod(trimString(tokens[0]));
            instance.translation[1] = std::stod(trimString(tokens[1]));
            instance.translation[2] = std::stod(trimString(tokens[2]));
        }
        lastPos = file.tellg();
    }

    m_instances.push_back(instance);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseEndInstance(std::ifstream& file)
{
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseEndAssembly(std::ifstream& file)
{
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::processData(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();

    // RgLog("Converting Abaqus data to FEModel...\n");

    // 计算总节点数和总单元数
    int totalNodes = 0;
    int totalElements = 0;
    for (const auto& instance : m_instances)
    {
        for (const auto& part : m_parts)
        {
            if (part.name == instance.partName)
            {
                totalNodes += part.nodes.size();
                totalElements += part.elements.size();
                break;
            }
        }
    }

    // 创建所有节点
    mesh.CreateNodes(totalNodes);

    // 初始化全局偏移
    m_globalNodeOffset = 0;
    m_globalElemOffset = 0;

    // Process each instance
    for (const auto& instance : m_instances)
    {
        // Find the corresponding part
        AbaqusPart* part = nullptr;
        for (auto& p : m_parts)
        {
            if (p.name == instance.partName)
            {
                part = &p;
                break;
            }
        }

        if (!part)
        {
            m_lastError = "Part not found for instance: " + instance.partName;
            return false;
        }

        // Create domain for this part
        RgDomain* domain = createDomain(*part, &mesh);
        if (!domain)
            return false;

        // Create nodes with global mapping
        std::map<int, int> nodeMap;  // Abaqus局部ID -> FEM全局索引
        if (!createNodes(*part, domain, &mesh, nodeMap, instance))
            return false;

        // Create elements with global node mapping
        if (!createElements(*part, domain, nodeMap))
            return false;

        // Create node sets
        if (!createNodeSets(instance, &mesh, nodeMap))
            return false;

        // Create element sets
        if (!createElementSets(instance, &mesh, domain))
            return false;

        // Create surfaces
        if (!createSurfaces(*part, &mesh, domain))
            return false;

        // Add domain to mesh
        mesh.AddDomain(domain);

        // Update global offset for next instance
        m_globalNodeOffset += part->nodes.size();
        m_globalElemOffset += part->elements.size();
    }

    //// Create boundary conditions
    // if (!createBoundaryConditions(fem))
    //     return false;

    //// Create loads
    // if (!createLoads(fem))
    //     return false;

    // Update mesh bounding box
    mesh.UpdateBox();

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createNodes(const AbaqusPart& part, RgDomain* domain, FEMesh* mesh, std::map<int, int>& nodeMap,
                               const AbaqusInstance& instance)
{
    for (size_t i = 0; i < part.nodes.size(); i++)
    {
        const AbaqusNode& abqNode = part.nodes[i];

        // 计算全局节点索引
        int globalIndex = m_globalNodeOffset + i;

        // 建立Abaqus局部ID到FEM全局索引的映射
        nodeMap[abqNode.id] = globalIndex;

        // 获取mesh中的节点
        FENode& meshNode = mesh->Node(globalIndex);

        // 应用instance的变换（平移和旋转）
        Vector3d pos(abqNode.x, abqNode.y, abqNode.z);

        // 应用平移
        pos.x += instance.translation[0];
        pos.y += instance.translation[1];
        pos.z += instance.translation[2];

        // TODO: 如果需要，应用旋转变换
        // applyRotation(pos, instance.rotation);

        meshNode.m_r0 = pos;
        meshNode.m_rt = pos;

        meshNode.SetID(globalIndex + 1);  // ID从1开始

        // 添加到全局映射
        m_globalNodeMap[abqNode.id] = globalIndex;

        domain->AddNode(&meshNode);
    }

    m_totalNodes += part.nodes.size();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createElements(const AbaqusPart& part, RgDomain* domain, const std::map<int, int>& nodeMap)
{
    // 为domain分配单元
    domain->Create(part.elements.size(), ElementType::FE_HEX8G8);

    for (size_t i = 0; i < part.elements.size(); i++)
    {
        const AbaqusElement& abqElem = part.elements[i];

        // 根据单元类型创建不同的单元
        RgElement* elem = nullptr;

        switch (abqElem.type)
        {
            case FE_HEX8:  // 假设您已经定义了这些常量
            {
                RgHex8Element* hex8 = dynamic_cast<RgHex8Element*>(&domain->ElementRef(i));
                if (hex8)
                {
                    elem = hex8;

                    // 设置节点连接性（使用全局节点索引）
                    if (abqElem.nodes.size() != 8)
                    {
                        m_lastError = "HEX8 element must have 8 nodes";
                        return false;
                    }

                    for (int j = 0; j < 8; j++)
                    {
                        int abaqusNodeId = abqElem.nodes[j];
                        auto it = nodeMap.find(abaqusNodeId);
                        if (it == nodeMap.end())
                        {
                            m_lastError = "Node ID not found: " + std::to_string(abaqusNodeId);
                            return false;
                        }

                        // 设置单元的节点索引
                        hex8->setNodeId(j, it->second);
                    }
                }
            }
            break;

            case FE_TET4:
            {
                // 类似处理TET4单元
                // RgTet4Element* tet4 = ...
            }
            break;

                // 添加其他单元类型...

            default:
                // RgLogWarning("Unsupported element type: %d", abqElem.type);
                continue;
        }

        if (elem)
        {
            // 设置单元ID（全局编号）
            int globalElemId = m_globalElemOffset + i;
            elem->setId(globalElemId + 1);  // ID从1开始

            // 添加到全局映射
            m_globalElemMap[abqElem.id] = globalElemId;
        }
    }

    m_totalElements += part.elements.size();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createNodeSets(const AbaqusInstance& part, FEMesh* mesh, const std::map<int, int>& nodeMap)
{
    for (const auto& nset : part.nodeSets)
    {
        FENodeSet* nodeSet = new FENodeSet(mesh->GetFEModel());
        nodeSet->SetName(nset.first);

        std::vector<int> nodeIndices;
        for (int abqNodeId : nset.second)
        {
            auto it = nodeMap.find(abqNodeId);
            if (it != nodeMap.end())
            {
                // 使用全局节点索引
                nodeIndices.push_back(it->second);
            }
            else
            {
                // RgLogWarning("Node ID %d not found in node set %s", abqNodeId, nset.first.c_str());
            }
        }

        // 根据您的FENodeSet API设置节点
        // 假设有如下方法（需要根据实际API调整）:
        // nodeSet->SetNodes(nodeIndices);
        // 或者
        // for (int idx : nodeIndices) {
        //     nodeSet->AddNode(idx);
        // }

        mesh->AddNodeSet(nodeSet);
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createElementSets(const AbaqusInstance& part, FEMesh* mesh, RgDomain* domain)
{
    for (const auto& elset : part.elementSets)
    {
        RgElementSet* elemSet = new RgElementSet(mesh->GetFEModel());
        elemSet->SetName(elset.first);

        std::vector<int> elemIndices;
        for (int abqElemId : elset.second)
        {
            //// 在当前part的elements中查找
            // for (size_t i = 0; i < part.elements.size(); i++)
            //{
            //     if (part.elements[i].id == abqElemId)
            //     {
            //         int globalElemIndex = m_globalElemOffset - part.elements.size() + i;
            //         elemIndices.push_back(globalElemIndex);
            //         break;
            //     }
            // }
        }

        // 根据您的RgElementSet API设置单元
        // elemSet->SetElements(elemIndices);
        // 或者从domain中添加
        // for (int idx : elemIndices) {
        //     elemSet->AddElement(&domain->ElementRef(idx));
        // }

        mesh->AddElementSet(elemSet);
    }

    return true;
}

//-----------------------------------------------------------------------------
RgDomain* AbaqusImport::createDomain(const AbaqusPart& part, FEMesh* mesh)
{
    // 根据part中的单元类型创建对应的domain
    // 这里假设一个part中的所有单元类型相同

    if (part.elements.empty())
    {
        m_lastError = "Part has no elements";
        return nullptr;
    }

    // 检查第一个单元的类型来确定domain类型
    int elemType = part.elements[0].type;

    RgDomain* domain = nullptr;

    // 根据单元类型创建domain
    if (elemType == FE_HEX8 || elemType == FE_HEX20 || elemType == FE_HEX27 || elemType == FE_TET4 ||
        elemType == FE_TET10 || elemType == FE_PENTA6 || elemType == FE_PENTA15)
    {
        // 实体单元 -> 实体域
        domain = new RgSolidDomain(mesh->GetFEModel());
        // RgLog("  Creating solid domain: %s\n", part.name.c_str());
    }
    else if (elemType == FE_QUAD4 || elemType == FE_QUAD8 || elemType == FE_TRI3 || elemType == FE_TRI6)
    {
        // 壳单元 -> 壳域
        // domain = new RgShellDomain(mesh->GetFEModel());
        // RgLog("  Creating shell domain: %s\n", part.name.c_str());
    }
    else if (elemType == FE_LINE2)
    {
        // 梁/桁架单元 -> 桁架域
        // domain = new RgTrussDomain(mesh->GetFEModel());
        // RgLog("  Creating truss domain: %s\n", part.name.c_str());
    }
    else
    {
        m_lastError = "Unsupported element type for domain creation";
        return nullptr;
    }

    if (domain)
    {
        domain->SetName(part.name);
    }

    return domain;
}

//-----------------------------------------------------------------------------
// 辅助函数：应用旋转变换（如果需要）
void AbaqusImport::applyRotation(Vector3d& pos, const double rotation[7])
{
    // 如果rotation定义为轴-角表示：rotation[0-2]是轴，rotation[3]是角度
    // 或者rotation[0-6]是四元数
    // 根据Abaqus的具体格式实现旋转

    // 示例：如果是轴-角表示
    // Vector3d axis(rotation[0], rotation[1], rotation[2]);
    // double angle = rotation[3];
    // 应用旋转矩阵...

    // TODO: 根据实际需要实现
}

// ============================================================================
// convertElementType函数需要返回实际的枚举值
// ============================================================================


int AbaqusImport::convertElementType(const std::string& abqType)
{
    std::string type = toUpper(abqType);

    // Solid elements
    if (type == "C3D8" || type == "C3D8R" || type == "C3D8I")
        return FE_HEX8;
    if (type == "C3D20" || type == "C3D20R")
        return FE_HEX20;
    if (type == "C3D27")
        return FE_HEX27;
    if (type == "C3D4")
        return FE_TET4;
    if (type == "C3D10" || type == "C3D10M")
        return FE_TET10;
    if (type == "C3D15")
        return FE_PENTA15;
    if (type == "C3D6")
        return FE_PENTA6;

    // Shell elements
    if (type == "S4" || type == "S4R")
        return FE_QUAD4;
    if (type == "S8R" || type == "S8R5")
        return FE_QUAD8;
    if (type == "S3" || type == "S3R")
        return FE_TRI3;
    if (type == "S6")
        return FE_TRI6;

    // Beam/Truss
    if (type == "T3D2" || type == "T2D2")
        return FE_LINE2;
    if (type == "B31" || type == "B32")
        return FE_LINE2;

    // RgLogWarning("Unknown element type: %s", abqType.c_str());
    return -1;
}

// ============================================================================
// getElementNodeCount函数需要返回实际的节点数
// ============================================================================

int AbaqusImport::getElementNodeCount(int elemType)
{
    switch (elemType)
    {
        case FE_HEX8:
            return 8;
        case FE_HEX20:
            return 20;
        case FE_HEX27:
            return 27;
        case FE_TET4:
            return 4;
        case FE_TET10:
            return 10;
        case FE_PENTA6:
            return 6;
        case FE_PENTA15:
            return 15;
        case FE_QUAD4:
            return 4;
        case FE_QUAD8:
            return 8;
        case FE_TRI3:
            return 3;
        case FE_TRI6:
            return 6;
        case FE_LINE2:
            return 2;
        default:
            return 0;
    }
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createSurfaces(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain)
{
    for (const auto& surf : part.surfaces)
    {
        FEFacetSet* facetSet = new FEFacetSet(mesh->GetFEModel());
        facetSet->SetName(surf.first);

        // TODO: Create facets from surface data
        mesh->AddFacetSet(facetSet);
    }

    return true;
}

////=============================================================================
//// Create Boundary Conditions
////=============================================================================
//
// bool AbaqusImport::createBoundaryConditions(FEModel* fem)
//{
//    RgLog("\n");
//    RgLog("=================================================\n");
//    RgLog("Creating Boundary Conditions\n");
//    RgLog("=================================================\n");
//
//    if (m_boundaryData.empty())
//    {
//        RgLogWarning("No boundary conditions to create\n");
//        return true;  // Not an error, just no BCs
//    }
//
//    FEMesh& mesh = fem->GetMesh();
//    int successCount = 0;
//    int failCount = 0;
//
//    // Process each parsed boundary condition
//    for (size_t i = 0; i < m_boundaryData.size(); ++i)
//    {
//        const BoundaryData& bcData = m_boundaryData[i];
//
//        // Find the node set
//        FENodeSet* nodeSet = mesh.FindNodeSet(bcData.nodeSetName);
//        if (!nodeSet)
//        {
//            RgLogError("Node set '%s' not found for boundary condition %d", bcData.nodeSetName.c_str(), (int)i + 1);
//            failCount++;
//            continue;
//        }
//
//        // Create appropriate boundary condition based on type
//        RgBoundaryCondition* bc = nullptr;
//
//        if (bcData.isEncastre)
//        {
//            // ENCASTRE - fix all DOFs
//            bc = RgBCFactory::CreateEncastre(fem, nodeSet);
//            bc->SetName(bcData.name.empty() ? bcData.nodeSetName + "_ENCASTRE" : bcData.name);
//
//            RgLog("  [%d] Created ENCASTRE BC on node set '%s'\n", successCount + 1, bcData.nodeSetName.c_str());
//        }
//        else
//        {
//            // Check if it's a fixed BC or prescribed BC
//            bool isFixed = (fabs(bcData.magnitude) < 1e-20);
//
//            // Create BC for each DOF in the range
//            for (int dof = bcData.firstDOF; dof <= bcData.lastDOF; ++dof)
//            {
//                int internalDOF = ConvertAbaqusDOF(dof);
//
//                if (isFixed)
//                {
//                    // Fixed BC (zero displacement/rotation)
//                    int dofMask = (1 << internalDOF);
//                    RgFixedBC* fixedBC = RgBCFactory::CreateFixed(fem, nodeSet, dofMask);
//
//                    std::string name = bcData.name.empty() ? bcData.nodeSetName + "_Fixed_DOF" + std::to_string(dof)
//                                                           : bcData.name + "_DOF" + std::to_string(dof);
//                    fixedBC->SetName(name);
//
//                    fem->AddBoundaryCondition(fixedBC);
//
//                    RgLog("  [%d] Created Fixed BC on '%s', DOF %d (internal %d)\n", successCount + 1,
//                          bcData.nodeSetName.c_str(), dof, internalDOF);
//                }
//                else
//                {
//                    // Prescribed BC (non-zero value)
//                    if (internalDOF < 3)
//                    {
//                    // Displacement
//                    RgPrescribedDisplacement* presBC =
//                        RgBCFactory::CreatePrescribed(fem, nodeSet, internalDOF, bcData.magnitude);
//
//                    std::string name = bcData.name.empty() ? bcData.nodeSetName + "_Disp_DOF" + std::to_string(dof)
//                                                           : bcData.name + "_DOF" + std::to_string(dof);
//                    presBC->SetName(name);
//
//                    fem->AddBoundaryCondition(presBC);
//
//                    RgLog("  [%d] Created Prescribed Displacement BC on '%s', "
//                          "DOF %d = %g\n",
//                          successCount + 1, bcData.nodeSetName.c_str(), dof, bcData.magnitude);
//                    }
//                    else
//                    {
//                    // Rotation
//                    RgPrescribedRotation* rotBC = new RgPrescribedRotation(fem);
//                    rotBC->SetNodeSet(nodeSet);
//                    rotBC->SetDOF(internalDOF);
//                    rotBC->SetScale(bcData.magnitude);
//
//                    std::string name = bcData.name.empty() ? bcData.nodeSetName + "_Rot_DOF" + std::to_string(dof)
//                                                           : bcData.name + "_DOF" + std::to_string(dof);
//                    rotBC->SetName(name);
//
//                    fem->AddBoundaryCondition(rotBC);
//
//                    RgLog("  [%d] Created Prescribed Rotation BC on '%s', "
//                          "DOF %d = %g\n",
//                          successCount + 1, bcData.nodeSetName.c_str(), dof, bcData.magnitude);
//                    }
//                }
//
//                successCount++;
//            }
//
//            continue;  // Skip the successCount++ at the end
//        }
//
//        // Add the BC to the model
//        if (bc)
//        {
//            fem->AddBoundaryCondition(bc);
//            successCount++;
//        }
//    }
//
//    // Summary
//    RgLog("-------------------------------------------------\n");
//    RgLog("Boundary Conditions Summary:\n");
//    RgLog("  Total entries parsed: %d\n", (int)m_boundaryData.size());
//    RgLog("  Successfully created: %d\n", successCount);
//    RgLog("  Failed: %d\n", failCount);
//    RgLog("  Total BCs in model: %d\n", fem->BoundaryConditions());
//    RgLog("=================================================\n\n");
//
//    return (failCount == 0);
//}
//
////=============================================================================
//// Create Loads
////=============================================================================
//
// bool AbaqusImport::createLoads(FEModel* fem)
//{
//    RgLog("\n");
//    RgLog("=================================================\n");
//    RgLog("Creating Loads\n");
//    RgLog("=================================================\n");
//
//    if (m_loadData.empty())
//    {
//        RgLogWarning("No loads to create\n");
//        return true;  // Not an error, just no loads
//    }
//
//    FEMesh& mesh = fem->GetMesh();
//    int successCount = 0;
//    int failCount = 0;
//
//    // Process each parsed load
//    for (size_t i = 0; i < m_loadData.size(); ++i)
//    {
//        const LoadData& loadData = m_loadData[i];
//
//        RgLoad* load = nullptr;
//
//        if (loadData.isConcentrated)
//        {
//            // Concentrated (nodal) load
//            load = CreateConcentratedLoad(fem, loadData);
//            if (load)
//            {
//                RgLog("  [%d] Created concentrated load on '%s'\n", successCount + 1, loadData.targetName.c_str());
//            }
//        }
//        else if (loadData.isDistributed)
//        {
//            // Distributed load
//            load = CreateDistributedLoad(fem, loadData);
//            if (load)
//            {
//                RgLog("  [%d] Created distributed load '%s' on '%s'\n", successCount + 1, loadData.loadType.c_str(),
//                      loadData.targetName.c_str());
//            }
//        }
//
//        if (load)
//        {
//            // Set name
//            std::string name = loadData.name.empty() ? loadData.targetName + "_Load" : loadData.name;
//            load->SetName(name);
//
//            // Add to model
//            fem->AddModelLoad(load);
//            successCount++;
//        }
//        else
//        {
//            RgLogError("Failed to create load %d", (int)i + 1);
//            failCount++;
//        }
//    }
//
//    // Summary
//    RgLog("-------------------------------------------------\n");
//    RgLog("Loads Summary:\n");
//    RgLog("  Total entries parsed: %d\n", (int)m_loadData.size());
//    RgLog("  Successfully created: %d\n", successCount);
//    RgLog("  Failed: %d\n", failCount);
//    RgLog("  Total loads in model: %d\n", fem->ModelLoads());
//    RgLog("=================================================\n\n");
//
//    return (failCount == 0);
//}
//
////=============================================================================
//// Helper: Create Concentrated Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreateConcentratedLoad(FEModel* fem, const LoadData& data)
//{
//    FEMesh& mesh = fem->GetMesh();
//
//    // Find node set
//    FENodeSet* nodeSet = mesh.FindNodeSet(data.targetName);
//    if (!nodeSet)
//    {
//        RgLogError("Node set '%s' not found for concentrated load", data.targetName.c_str());
//        return nullptr;
//    }
//
//    int internalDOF = ConvertAbaqusDOF(data.dof);
//
//    if (internalDOF < 3)
//    {
//        // Force
//        RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(fem, nodeSet, internalDOF, data.magnitude);
//
//        RgLog("    Type: Nodal Force, DOF %d (Abaqus DOF %d), Magnitude %g\n", internalDOF, data.dof, data.magnitude);
//
//        return load;
//    }
//    else
//    {
//        // Moment
//        RgMomentLoad* load = new RgMomentLoad(fem);
//        load->SetNodeSet(nodeSet);
//
//        Vector3d moment(0, 0, 0);
//        moment[internalDOF - 3] = data.magnitude;
//        load->SetMoment(moment);
//        load->SetMagnitude(1.0);
//
//        RgLog("    Type: Moment, DOF %d (Abaqus DOF %d), Magnitude %g\n", internalDOF, data.dof, data.magnitude);
//
//        return load;
//    }
//}
//
////=============================================================================
//// Helper: Create Distributed Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreateDistributedLoad(FEModel* fem, const LoadData& data)
//{
//    FEMesh& mesh = fem->GetMesh();
//    std::string loadType = data.loadType;
//
//    // Convert to uppercase for comparison
//    std::transform(loadType.begin(), loadType.end(), loadType.begin(), ::toupper);
//
//    if (loadType == "P")
//    {
//        // Pressure load
//        return CreatePressureLoad(fem, data);
//    }
//    else if (loadType == "TRVEC")
//    {
//        // Traction vector
//        return CreateTractionLoad(fem, data);
//    }
//    else if (loadType == "GRAV")
//    {
//        // Gravity
//        return CreateGravityLoad(fem, data);
//    }
//    else if (loadType == "CENTRIF")
//    {
//        // Centrifugal
//        return CreateCentrifugalLoad(fem, data);
//    }
//    else
//    {
//        RgLogError("Unknown distributed load type: %s", loadType.c_str());
//        return nullptr;
//    }
//}
//
////=============================================================================
//// Helper: Create Pressure Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreatePressureLoad(FEModel* fem, const LoadData& data)
//{
//    FEMesh& mesh = fem->GetMesh();
//
//    // Try to find as facet set first
//    FEFacetSet* facetSet = mesh.FindFacetSet(data.targetName);
//    if (facetSet)
//    {
//        RgSurfaceLoad* load = RgLoadFactory::CreatePressure(fem, facetSet, data.magnitude);
//
//        RgLog("    Type: Pressure, Magnitude %g Pa\n", data.magnitude);
//
//        return load;
//    }
//
//    // Try to find as surface
//    FESurface* surface = mesh.FindSurface(data.targetName);
//    if (surface)
//    {
//        RgSurfaceLoad* load = RgLoadFactory::CreatePressure(fem, surface, data.magnitude);
//
//        RgLog("    Type: Pressure, Magnitude %g Pa\n", data.magnitude);
//
//        return load;
//    }
//
//    RgLogError("Surface/Facet set '%s' not found for pressure load", data.targetName.c_str());
//    return nullptr;
//}
//
////=============================================================================
//// Helper: Create Traction Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreateTractionLoad(FEModel* fem, const LoadData& data)
//{
//    FEMesh& mesh = fem->GetMesh();
//
//    // Find facet set
//    FEFacetSet* facetSet = mesh.FindFacetSet(data.targetName);
//    if (!facetSet)
//    {
//        RgLogError("Facet set '%s' not found for traction load", data.targetName.c_str());
//        return nullptr;
//    }
//
//    // Create surface from facet set
//    FESurface* surface = mesh.CreateSurface(*facetSet);
//    if (!surface)
//    {
//        RgLogError("Failed to create surface from facet set '%s'", data.targetName.c_str());
//        return nullptr;
//    }
//
//    Vector3d traction(data.x, data.y, data.z);
//
//    RgSurfaceLoad* load = RgLoadFactory::CreateTraction(fem, surface, traction);
//    load->SetMagnitude(data.magnitude);
//
//    RgLog("    Type: Traction, Direction (%g, %g, %g), Magnitude %g\n", data.x, data.y, data.z, data.magnitude);
//
//    return load;
//}
//
////=============================================================================
//// Helper: Create Gravity Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreateGravityLoad(FEModel* fem, const LoadData& data)
//{
//    // Create gravity vector
//    Vector3d direction(data.x, data.y, data.z);
//    direction.unit();  // Normalize
//
//    Vector3d g = direction * data.magnitude;
//
//    RgBodyLoad* load = RgLoadFactory::CreateGravity(fem, g);
//
//    RgLog("    Type: Gravity, g = (%g, %g, %g) m/s²\n", g.x, g.y, g.z);
//
//    return load;
//}
//
////=============================================================================
//// Helper: Create Centrifugal Load
////=============================================================================
//
// RgLoad* AbaqusImport::CreateCentrifugalLoad(FEModel* fem, const LoadData& data)
//{
//    // In Abaqus, CENTRIF format is:
//    // name, CENTRIF, omega^2, x0, y0, z0, x1, y1, z1
//    // where (x0,y0,z0) is a point on the axis
//    // and (x1,y1,z1) is the axis direction
//
//    Vector3d origin(data.x, data.y, data.z);
//
//    // For now, we'll use a default Z-axis
//    // In a complete implementation, you'd need to extend LoadData
//    // to store the axis direction separately
//    Vector3d axis(0, 0, 1);
//
//    // TODO: Parse axis direction from additional data
//    // This would require extending the LoadData structure
//
//    // Abaqus stores omega^2, we need omega
//    double omega = sqrt(fabs(data.magnitude));
//
//    RgBodyLoad* load = RgLoadFactory::CreateCentrifugal(fem, axis, origin, omega);
//
//    RgLog("    Type: Centrifugal, omega = %g rad/s, origin = (%g, %g, %g)\n", omega, origin.x, origin.y, origin.z);
//    RgLogWarning("    Note: Axis direction defaulted to Z-axis (0,0,1)\n");
//
//    return load;
//}
//
////=============================================================================
//// Helper: Convert Abaqus DOF to Internal DOF
////=============================================================================
//
// int AbaqusImport::ConvertAbaqusDOF(int abaqusDOF)
//{
//    // Abaqus: 1,2,3 = X,Y,Z displacement; 4,5,6 = X,Y,Z rotation
//    // Internal: 0,1,2 = X,Y,Z displacement; 3,4,5 = X,Y,Z rotation
//
//    if (abaqusDOF < 1 || abaqusDOF > 6)
//    {
//        RgLogWarning("Invalid Abaqus DOF: %d, using 0", abaqusDOF);
//        return 0;
//    }
//
//    return abaqusDOF - 1;
//}
//
////=============================================================================
//// Optional: Attach Load Controllers
////=============================================================================
//
// bool AbaqusImport::attachLoadControllers(FEModel* fem)
//{
//    RgLog("\n");
//    RgLog("=================================================\n");
//    RgLog("Attaching Load Controllers\n");
//    RgLog("=================================================\n");
//
//    // This is optional - you can attach load curves to BCs and loads
//    // after they've been created
//
//    // Example: Create a default ramp load curve
//    if (fem->LoadControllers() == 0)
//    {
//        RgLoadCurve* defaultCurve = new RgLoadCurve();
//        defaultCurve->AddPoint(0.0, 0.0);
//        defaultCurve->AddPoint(1.0, 1.0);
//        defaultCurve->SetName("DefaultRamp");
//        defaultCurve->SetInterpolation(RgLoadController::INTERP_LINEAR);
//        fem->AddLoadController(defaultCurve);
//
//        RgLog("  Created default ramp load curve\n");
//    }
//
//    // Attach to prescribed BCs if needed
//    int attachCount = 0;
//    int nbc = fem->BoundaryConditions();
//    for (int i = 0; i < nbc; ++i)
//    {
//        RgBoundaryCondition* bc = dynamic_cast<RgBoundaryCondition*>(fem->BoundaryCondition(i));
//
//        // Check if it's a prescribed BC with non-zero magnitude
//        RgPrescribedDisplacement* presBC = dynamic_cast<RgPrescribedDisplacement*>(bc);
//
//        if (presBC && fabs(presBC->GetScale()) > 1e-12)
//        {
//            // Attach default load curve if no controller set
//            if (!presBC->GetLoadController())
//            {
//                presBC->SetLoadController(dynamic_cast<RgLoadController*>(fem->GetLoadController(0)));
//                attachCount++;
//            }
//        }
//    }
//
//    if (attachCount > 0)
//    {
//        RgLog("  Attached load controller to %d prescribed BCs\n", attachCount);
//    }
//
//    // Attach to loads if needed
//    attachCount = 0;
//    int nloads = fem->ModelLoads();
//    for (int i = 0; i < nloads; ++i)
//    {
//        RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
//
//        if (load && !load->GetLoadController())
//        {
//            load->SetLoadController(dynamic_cast<RgLoadController*>(fem->GetLoadController(0)));
//            attachCount++;
//        }
//    }
//
//    if (attachCount > 0)
//    {
//        RgLog("  Attached load controller to %d loads\n", attachCount);
//    }
//
//    RgLog("=================================================\n\n");
//
//    return true;
//}
//
////=============================================================================
//// Validation
////=============================================================================
//
// bool AbaqusImport::validateBCsAndLoads(FEModel* fem)
//{
//    RgLog("\n");
//    RgLog("=================================================\n");
//    RgLog("Validating Boundary Conditions and Loads\n");
//    RgLog("=================================================\n");
//
//    bool success = true;
//    FEMesh& mesh = fem->GetMesh();
//
//    // Validate BCs
//    int nbc = fem->BoundaryConditions();
//    RgLog("Checking %d boundary conditions...\n", nbc);
//
//    for (int i = 0; i < nbc; ++i)
//    {
//        RgBoundaryCondition* bc = dynamic_cast<RgBoundaryCondition*>(fem->BoundaryCondition(i));
//
//        if (!bc)
//        {
//            RgLogError("  BC %d: Invalid BC pointer", i);
//            success = false;
//            continue;
//        }
//
//        FENodeSet* nodeSet = bc->GetNodeSet();
//        if (!nodeSet)
//        {
//            RgLogError("  BC %d (%s): No node set assigned", i, bc->GetName().c_str());
//            success = false;
//            continue;
//        }
//
//        if (nodeSet->Size() == 0)
//        {
//            RgLogWarning("  BC %d (%s): Node set '%s' is empty", i, bc->GetName().c_str(),
//            nodeSet->GetName().c_str());
//        }
//    }
//
//    // Validate Loads
//    int nloads = fem->ModelLoads();
//    RgLog("Checking %d loads...\n", nloads);
//
//    for (int i = 0; i < nloads; ++i)
//    {
//        RgLoad* load = dynamic_cast<RgLoad*>(fem->ModelLoad(i));
//
//        if (!load)
//        {
//            RgLogError("  Load %d: Invalid load pointer", i);
//            success = false;
//            continue;
//        }
//
//        // Check nodal loads
//        RgNodalLoad* nodalLoad = dynamic_cast<RgNodalLoad*>(load);
//        if (nodalLoad)
//        {
//            FENodeSet* nodeSet = nodalLoad->GetNodeSet();
//            if (!nodeSet)
//            {
//                RgLogError("  Load %d (%s): No node set assigned", i, load->GetName().c_str());
//                success = false;
//            }
//            else if (nodeSet->Size() == 0)
//            {
//                RgLogWarning("  Load %d (%s): Node set '%s' is empty", i, load->GetName().c_str(),
//                             nodeSet->GetName().c_str());
//            }
//        }
//
//        // Check surface loads
//        RgSurfaceLoad* surfLoad = dynamic_cast<RgSurfaceLoad*>(load);
//        if (surfLoad)
//        {
//            if (!surfLoad->GetSurface() && !surfLoad->GetFacetSet())
//            {
//                RgLogError("  Load %d (%s): No surface/facet set assigned", i, load->GetName().c_str());
//                success = false;
//            }
//        }
//    }
//
//    if (success)
//    {
//        RgLog("  All boundary conditions and loads are valid\n");
//    }
//    else
//    {
//        RgLogError("  Validation failed - some BCs/loads are invalid\n");
//    }
//
//    RgLog("=================================================\n\n");
//
//    return success;
//}


//-----------------------------------------------------------------------------
// int AbaqusImport::convertElementType(const std::string& abqType)
//{
//    std::string type = toUpper(abqType);
//
//    //// Solid elements
//    //if (type == "C3D8" || type == "C3D8R" || type == "C3D8I")
//    //    return FE_HEX8;
//    //if (type == "C3D20" || type == "C3D20R")
//    //    return FE_HEX20;
//    //if (type == "C3D27")
//    //    return FE_HEX27;
//    //if (type == "C3D4")
//    //    return FE_TET4;
//    //if (type == "C3D10" || type == "C3D10M")
//    //    return FE_TET10;
//    //if (type == "C3D15")
//    //    return FE_PENTA15;
//    //if (type == "C3D6")
//    //    return FE_PENTA6;
//
//    //// Shell elements
//    //if (type == "S4" || type == "S4R")
//    //    return FE_QUAD4;
//    //if (type == "S8R" || type == "S8R5")
//    //    return FE_QUAD8;
//    //if (type == "S3" || type == "S3R")
//    //    return FE_TRI3;
//    //if (type == "S6")
//    //    return FE_TRI6;
//
//    //// Beam/Truss
//    //if (type == "T3D2" || type == "T2D2")
//    //    return FE_LINE2;
//    //if (type == "B31" || type == "B32")
//    //    return FE_LINE2;
//
//    ////RgLogWarning("Unknown element type: %s", abqType.c_str());
//    return -1;
//}

////-----------------------------------------------------------------------------
// int AbaqusImport::getElementNodeCount(int elemType)
//{
//     switch (elemType)
//     {
//         /*case FE_HEX8:
//             return 8;
//         case FE_HEX20:
//             return 20;
//         case FE_HEX27:
//             return 27;
//         case FE_TET4:
//             return 4;
//         case FE_TET10:
//             return 10;
//         case FE_PENTA6:
//             return 6;
//         case FE_PENTA15:
//             return 15;
//         case FE_QUAD4:
//             return 4;
//         case FE_QUAD8:
//             return 8;
//         case FE_TRI3:
//             return 3;
//         case FE_TRI6:
//             return 6;
//         case FE_LINE2:
//             return 2;*/
//     default:
//         return 0;
//     }
// }

//-----------------------------------------------------------------------------
bool AbaqusImport::validatePart(const AbaqusPart& part)
{
    if (part.name.empty())
    {
        m_lastError = "Part must have a name";
        return false;
    }

    if (part.nodes.empty())
    {
        m_lastError = "Part has no nodes";
        return false;
    }

    if (part.elements.empty())
    {
        m_lastError = "Part has no elements";
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::validateNodeId(int nodeId, const std::map<int, int>& nodeMap)
{
    return (nodeMap.find(nodeId) != nodeMap.end());
}

//-----------------------------------------------------------------------------
std::string AbaqusImport::readKeywordLine(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);
    m_lineNumber++;
    return line;
}

//-----------------------------------------------------------------------------
void AbaqusImport::parseKeywordParams(const std::string& line, std::map<std::string, std::string>& params)
{
    std::vector<std::string> tokens = splitString(line, ',');

    for (size_t i = 1; i < tokens.size(); i++)
    {
        std::string token = trimString(tokens[i]);
        size_t pos = token.find('=');

        if (pos != std::string::npos)
        {
            std::string key = toUpper(trimString(token.substr(0, pos)));
            std::string value = trimString(token.substr(pos + 1));
            params[key] = value;
        }
    }
}

//-----------------------------------------------------------------------------
std::vector<std::string> AbaqusImport::splitString(const std::string& str, char delimiter)
{
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter))
    {
        tokens.push_back(token);
    }

    return tokens;
}

//-----------------------------------------------------------------------------
std::string AbaqusImport::trimString(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

//-----------------------------------------------------------------------------
std::string AbaqusImport::toUpper(const std::string& str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::isKeyword(const std::string& line)
{
    if (line.empty())
        return false;
    return (line[0] == '*' && (line.size() < 2 || line[1] != '*'));
}

//-----------------------------------------------------------------------------
bool AbaqusImport::skipToNextKeyword(std::ifstream& file)
{
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);
        if (!line.empty() && line[0] == '*')
        {
            file.seekg(lastPos);
            m_lineNumber--;
            return true;
        }
        lastPos = file.tellg();
    }
    return false;
}