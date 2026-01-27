#include "AbaqusImport.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENode.h"
#include "femcore/Domain/RgDomain.h"
#include "femcore/Domain/RgSolidDomain.h"
#include "elements/RgElementSet.h"
#include "femcore/FENodeSet.h"
#include "femcore/FEFacetSet.h"
#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEModelLoad.h"
#include "logger/log.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <iostream>

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
        //feLogError(m_lastError.c_str());
        return false;
    }

    //feLog("Reading Abaqus INP file: %s\n", filename);

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
        //feLogError("Failed to parse INP file: %s", m_lastError.c_str());
        return false;
    }

    // Convert parsed data to FEModel
    if (!processData(fem))
    {
        //feLogError("Failed to process INP data: %s", m_lastError.c_str());
        return false;
    }

    //feLog("Successfully imported %d nodes and %d elements from %d parts\n", m_totalNodes, m_totalElements,
    //      (int)m_parts.size());

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
                m_inAssembly = false;
            }
            else if (keyword.find("MATERIAL") == 0)
            {
                parseMaterial(file, line);
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
    //feLog("  Parsing part: %s\n", part.name.c_str());

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

                //feLog("    Section: ELSET=%s, MATERIAL=%s\n", elsetName.c_str(), materialName.c_str());

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

    //feLog("    Read %d nodes\n", nodeCount);
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

    //feLog("    Read %d %s elements\n", elemCount, elemType.c_str());
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
    //feLog("    Read node set '%s' with %d nodes\n", nsetName.c_str(), (int)nodeIds.size());
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
    //feLog("    Read element set '%s' with %d elements\n", elsetName.c_str(), (int)elemIds.size());
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
            faceData.push_back(0);  // Placeholder for element/elset reference
            faceData.push_back(std::stoi(trimString(tokens[1])));  // Face ID
            surfaceData.push_back(faceData);
        }
        lastPos = file.tellg();
    }

    part.surfaces[surfaceName] = surfaceData;
    //feLog("    Read surface '%s' with %d faces\n", surfaceName.c_str(), (int)surfaceData.size());
    return true;
}

//-----------------------------------------------------------------------------
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
    //feLog("  Parsing material: %s\n", material.name.c_str());

    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
        {
            lastPos = file.tellg();
            continue;
        }

        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            if (keyword.find("ELASTIC") == 0)
            {
                // Read elastic properties
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
                            material.type = "ELASTIC";
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
                        densityProp.push_back(std::stod(trimString(line)));
                        material.properties["DENSITY"] = densityProp;
                        lastPos = file.tellg();
                    }
                    else
                    {
                        file.seekg(lastPos);
                        m_lineNumber--;
                    }
                }
            }
            else
            {
                // Unknown material property, return to process next material/keyword
                file.seekg(lastPos);
                m_lineNumber--;
                break;
            }
        }
    }

    m_materials.push_back(material);
    m_currentMaterial.clear();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseStep(std::ifstream& file, const std::string& keywordLine)
{
    m_currentStep++;
    m_inStep = true;

    StepInfo step;
    step.name = "Step-" + std::to_string(m_currentStep + 1);

    //feLog("  Parsing step %d\n", m_currentStep + 1);

    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);

        if (line.empty() || (line.size() >= 2 && line[0] == '*' && line[1] == '*'))
        {
            lastPos = file.tellg();
            continue;
        }

        if (line[0] == '*')
        {
            std::string keyword = toUpper(line.substr(1));

            if (keyword.find("STATIC") == 0 || keyword.find("DYNAMIC") == 0)
            {
                step.procedure = (keyword.find("STATIC") == 0) ? "STATIC" : "DYNAMIC";
                // Read step parameters
                lastPos = file.tellg();
                if (std::getline(file, line))
                {
                    m_lineNumber++;
                    line = trimString(line);
                    if (!line.empty() && line[0] != '*')
                    {
                        std::vector<std::string> tokens = splitString(line, ',');
                        if (tokens.size() >= 1)
                            step.initialTimeIncrement = std::stod(trimString(tokens[0]));
                        if (tokens.size() >= 2)
                            step.timePeriod = std::stod(trimString(tokens[1]));
                        if (tokens.size() >= 3)
                            step.minTimeIncrement = std::stod(trimString(tokens[2]));
                        if (tokens.size() >= 4)
                            step.maxTimeIncrement = std::stod(trimString(tokens[3]));
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
                parseBoundary(file);
            }
            else if (keyword.find("CLOAD") == 0)
            {
                parseCload(file);
            }
            else if (keyword.find("DLOAD") == 0)
            {
                parseDload(file);
            }
            else if (keyword.find("END STEP") == 0)
            {
                break;
            }
            else
            {
                // Skip unknown keywords within step
                lastPos = file.tellg();
            }
        }
    }

    m_steps.push_back(step);
    m_inStep = false;
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseBoundary(std::ifstream& file)
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

        // Parse boundary condition: node_set, first_dof, last_dof, magnitude
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() >= 3)
        {
            BoundaryCondition bc;
            bc.nodeSetName = trimString(tokens[0]);
            bc.dof = std::stoi(trimString(tokens[1]));
            bc.value = (tokens.size() >= 4) ? std::stod(trimString(tokens[3])) : 0.0;
            bc.step = m_currentStep;
            m_boundaryConditions.push_back(bc);
        }
        lastPos = file.tellg();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseLoad(std::ifstream& file)
{
    // Generic load parsing - delegates to specific load types
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseCload(std::ifstream& file)
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

        // Parse concentrated load: node_set, dof, magnitude
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() >= 3)
        {
            ConcentratedLoad load;
            load.nodeSetName = trimString(tokens[0]);
            load.dof = std::stoi(trimString(tokens[1]));
            load.magnitude = std::stod(trimString(tokens[2]));
            load.step = m_currentStep;
            m_concentratedLoads.push_back(load);
        }
        lastPos = file.tellg();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::parseDload(std::ifstream& file)
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

        // Parse distributed load: surface, load_type, magnitude
        std::vector<std::string> tokens = splitString(line, ',');
        if (tokens.size() >= 3)
        {
            DistributedLoad load;
            load.surfaceName = trimString(tokens[0]);
            load.loadType = trimString(tokens[1]);
            load.magnitude = std::stod(trimString(tokens[2]));
            load.step = m_currentStep;
            m_distributedLoads.push_back(load);
        }
        lastPos = file.tellg();
    }

    return true;
}

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

        //feLog("Heading: %s\n", line.c_str());
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

    //feLog("  Parsing instance: %s (part: %s)\n", instance.name.c_str(), instance.partName.c_str());

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

    //feLog("Converting Abaqus data to FEModel...\n");

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
        if (!createNodes(*part, &mesh, nodeMap, instance))
            return false;

        // Create elements with global node mapping
        if (!createElements(*part, domain, nodeMap))
            return false;

        // Create node sets
        if (!createNodeSets(*part, &mesh, nodeMap))
            return false;

        // Create element sets
        if (!createElementSets(*part, &mesh, domain))
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

    // Create boundary conditions
    if (!createBoundaryConditions(fem))
        return false;

    // Create loads
    if (!createLoads(fem))
        return false;

    // Update mesh bounding box
    mesh.UpdateBox();

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createNodes(const AbaqusPart& part, FEMesh* mesh,
    std::map<int, int>& nodeMap,
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
    }

    m_totalNodes += part.nodes.size();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createElements(const AbaqusPart& part, RgDomain* domain,
    const std::map<int, int>& nodeMap)
{
    // 为domain分配单元
    domain->Create(part.elements.size());

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
            //feLogWarning("Unsupported element type: %d", abqElem.type);
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
bool AbaqusImport::createNodeSets(const AbaqusPart& part, FEMesh* mesh,
    const std::map<int, int>& nodeMap)
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
                //feLogWarning("Node ID %d not found in node set %s", abqNodeId, nset.first.c_str());
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
bool AbaqusImport::createElementSets(const AbaqusPart& part, FEMesh* mesh,
    RgDomain* domain)
{
    for (const auto& elset : part.elementSets)
    {
        RgElementSet* elemSet = new RgElementSet(mesh->GetFEModel());
        elemSet->SetName(elset.first);

        std::vector<int> elemIndices;
        for (int abqElemId : elset.second)
        {
            // 在当前part的elements中查找
            for (size_t i = 0; i < part.elements.size(); i++)
            {
                if (part.elements[i].id == abqElemId)
                {
                    int globalElemIndex = m_globalElemOffset - part.elements.size() + i;
                    elemIndices.push_back(globalElemIndex);
                    break;
                }
            }
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
    if (elemType == FE_HEX8 || elemType == FE_HEX20 || elemType == FE_HEX27 ||
        elemType == FE_TET4 || elemType == FE_TET10 ||
        elemType == FE_PENTA6 || elemType == FE_PENTA15)
    {
        // 实体单元 -> 实体域
        domain = new RgSolidDomain(mesh->GetFEModel());
        //feLog("  Creating solid domain: %s\n", part.name.c_str());
    }
    else if (elemType == FE_QUAD4 || elemType == FE_QUAD8 ||
        elemType == FE_TRI3 || elemType == FE_TRI6)
    {
        // 壳单元 -> 壳域
        // domain = new RgShellDomain(mesh->GetFEModel());
        //feLog("  Creating shell domain: %s\n", part.name.c_str());
    }
    else if (elemType == FE_LINE2)
    {
        // 梁/桁架单元 -> 桁架域
        // domain = new RgTrussDomain(mesh->GetFEModel());
        //feLog("  Creating truss domain: %s\n", part.name.c_str());
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
// 单元类型常量定义（如果还没有，需要在某个头文件中定义）
// ============================================================================


enum ElementType
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

    //feLogWarning("Unknown element type: %s", abqType.c_str());
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
bool AbaqusImport::processData(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();

    //feLog("Converting Abaqus data to FEModel...\n");

    // Process each part/instance
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

        // Create nodes
        std::map<int, int> nodeMap;
        if (!createNodes(*part, &mesh, nodeMap))
            return false;

        // Create elements
        if (!createElements(*part, domain, nodeMap))
            return false;

        // Create node sets
        if (!createNodeSets(*part, &mesh, nodeMap))
            return false;

        // Create element sets
        if (!createElementSets(*part, &mesh, domain))
            return false;

        // Create surfaces
        if (!createSurfaces(*part, &mesh, domain))
            return false;

        // Add domain to mesh
        mesh.AddDomain(domain);
    }

    // Create boundary conditions
    if (!createBoundaryConditions(fem))
        return false;

    // Create loads
    if (!createLoads(fem))
        return false;

    // Update mesh bounding box
    mesh.UpdateBox();

    return true;
}

//-----------------------------------------------------------------------------
RgDomain* AbaqusImport::createDomain(const AbaqusPart& part, FEMesh* mesh)
{
    // For now, create a solid domain
    // TODO: Detect domain type based on element types
    RgSolidDomain* domain = new RgSolidDomain(mesh->GetFEModel());
    domain->SetName(part.name);

    //feLog("  Creating domain: %s\n", part.name.c_str());
    return domain;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createNodes(const AbaqusPart& part, FEMesh* mesh, std::map<int, int>& nodeMap)
{
    int nodeOffset = m_totalNodes;

    for (const auto& node : part.nodes)
    {
        FENode feNode;
        feNode.m_r0 = Vector3d(node.x, node.y, node.z);
        feNode.m_rt = feNode.m_r0;
        feNode.SetID(node.id);

        // Map Abaqus node ID to FEModel node index
        nodeMap[node.id] = nodeOffset + (&node - &part.nodes[0]);
    }

    // Add nodes to mesh
    int currentNodeCount = mesh->Nodes();
    mesh->AddNodes(part.nodes.size());

    for (size_t i = 0; i < part.nodes.size(); i++)
    {
        FENode& meshNode = mesh->Node(currentNodeCount + i);
        meshNode.m_r0 = Vector3d(part.nodes[i].x, part.nodes[i].y, part.nodes[i].z);
        meshNode.m_rt = meshNode.m_r0;
        meshNode.SetID(part.nodes[i].id);
    }

    m_totalNodes += part.nodes.size();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createElements(const AbaqusPart& part, RgDomain* domain, const std::map<int, int>& nodeMap)
{
    // TODO: Create actual elements and add to domain
    // This requires knowledge of your element class hierarchy
    m_totalElements += part.elements.size();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createNodeSets(const AbaqusPart& part, FEMesh* mesh, const std::map<int, int>& nodeMap)
{
    for (const auto& nset : part.nodeSets)
    {
        FENodeSet* nodeSet = new FENodeSet(mesh->GetFEModel());
        nodeSet->SetName(nset.first);

        std::vector<int> nodeIndices;
        for (int abqId : nset.second)
        {
            auto it = nodeMap.find(abqId);
            if (it != nodeMap.end())
            {
                nodeIndices.push_back(it->second);
            }
        }

        // TODO: Set node indices to nodeSet based on your FENodeSet API
        mesh->AddNodeSet(nodeSet);
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createElementSets(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain)
{
    for (const auto& elset : part.elementSets)
    {
        RgElementSet* elemSet = new RgElementSet(mesh->GetFEModel());
        elemSet->SetName(elset.first);

        // TODO: Add elements to set based on your RgElementSet API
        mesh->AddElementSet(elemSet);
    }

    return true;
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

//-----------------------------------------------------------------------------
bool AbaqusImport::createBoundaryConditions(FEModel* fem)
{
    for (const auto& bc : m_boundaryConditions)
    {
        // TODO: Create FEBoundaryCondition objects
        // This depends on your boundary condition class implementation
        ////feLog("  Creating BC for node set '%s', DOF %d\n", bc.nodeSetName.c_str(), bc.dof);
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::createLoads(FEModel* fem)
{
    // Create concentrated loads
    for (const auto& load : m_concentratedLoads)
    {
        // TODO: Create concentrated load objects
        ////feLog("  Creating concentrated load for node set '%s', DOF %d, magnitude %g\n", load.nodeSetName.c_str(),
        //      load.dof, load.magnitude);
    }

    // Create distributed loads
    for (const auto& load : m_distributedLoads)
    {
        // TODO: Create distributed load (pressure) objects
        ////feLog("  Creating distributed load for surface '%s', type %s, magnitude %g\n", load.surfaceName.c_str(),
        //      load.loadType.c_str(), load.magnitude);
    }

    return true;
}

//-----------------------------------------------------------------------------
int AbaqusImport::convertElementType(const std::string& abqType)
{
    std::string type = toUpper(abqType);

    //// Solid elements
    //if (type == "C3D8" || type == "C3D8R" || type == "C3D8I")
    //    return FE_HEX8;
    //if (type == "C3D20" || type == "C3D20R")
    //    return FE_HEX20;
    //if (type == "C3D27")
    //    return FE_HEX27;
    //if (type == "C3D4")
    //    return FE_TET4;
    //if (type == "C3D10" || type == "C3D10M")
    //    return FE_TET10;
    //if (type == "C3D15")
    //    return FE_PENTA15;
    //if (type == "C3D6")
    //    return FE_PENTA6;

    //// Shell elements
    //if (type == "S4" || type == "S4R")
    //    return FE_QUAD4;
    //if (type == "S8R" || type == "S8R5")
    //    return FE_QUAD8;
    //if (type == "S3" || type == "S3R")
    //    return FE_TRI3;
    //if (type == "S6")
    //    return FE_TRI6;

    //// Beam/Truss
    //if (type == "T3D2" || type == "T2D2")
    //    return FE_LINE2;
    //if (type == "B31" || type == "B32")
    //    return FE_LINE2;

    ////feLogWarning("Unknown element type: %s", abqType.c_str());
    return -1;
}

//-----------------------------------------------------------------------------
int AbaqusImport::getElementNodeCount(int elemType)
{
    switch (elemType)
    {
        /*case FE_HEX8:
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
            return 2;*/
    default:
        return 0;
    }
}

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