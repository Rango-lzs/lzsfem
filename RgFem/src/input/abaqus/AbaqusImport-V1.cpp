#include "AbaqusImport.h"

#include "femcore/Domain/RgDomain.h"
#include "femcore/Domain/RgSolidDomain.h"
#include "elements/RgElementSet.h"
#include "femcore/FEFacetSet.h"
#include "femcore/FENodeSet.h"
#include "femcore/FEMesh.h"
#include "femcore/FEModel.h"
#include "femcore/FENode.h"
#include "logger/log.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>

//-----------------------------------------------------------------------------
AbaqusImport::AbaqusImport()
    : m_inAssembly(false)
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
bool AbaqusImport::Load(const char* filename, FEModel* fem)
{
    if (!filename || !fem)
    {
        m_lastError = "Invalid input parameters";
        feLogErrorEx(fem, m_lastError.c_str());
        return false;
    }

    feLogEx(fem, "Reading Abaqus INP file: %s\n", filename);

    // Clear previous data
    m_parts.clear();
    m_instances.clear();
    m_materials.clear();
    m_globalNodeMap.clear();
    m_globalElemMap.clear();
    m_totalNodes = 0;
    m_totalElements = 0;

    // Parse the INP file
    if (!ParseFile(filename))
    {
        feLogErrorEx(fem, "Failed to parse INP file: %s", m_lastError.c_str());
        return false;
    }

    // Convert parsed data to FEModel
    if (!ProcessData(fem))
    {
        feLogErrorEx(fem, "Failed to process INP data: %s", m_lastError.c_str());
        return false;
    }

    feLogEx(fem, "Successfully imported %d nodes and %d elements from %d parts\n", m_totalNodes, m_totalElements,
          (int)m_parts.size());

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseFile(const char* filename)
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
        line = TrimString(line);
        if (line.empty() || line[0] == '*' && line[1] == '*')
            continue;

        // Check if this is a keyword line
        if (line[0] == '*' && line[1] != '*')
        {
            std::string keyword = ToUpper(line.substr(1));

            // Parse based on keyword
            if (keyword.find("HEADING") == 0)
            {
                ParseHeading(file);
            }
            else if (keyword.find("PART") == 0)
            {
                if (!ParsePart(file))
                    return false;
            }
            else if (keyword.find("ASSEMBLY") == 0)
            {
                m_inAssembly = true;
                if (!ParseAssembly(file))
                    return false;
                m_inAssembly = false;
            }
            else if (keyword.find("MATERIAL") == 0)
            {
                ParseMaterial(file);
            }
            else if (keyword.find("STEP") == 0)
            {
                ParseStep(file);
            }
        }
    }

    file.close();
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParsePart(std::ifstream& file)
{
    std::string line = ReadKeywordLine(file);
    std::map<std::string, std::string> params;
    ParseKeywordParams(line, params);

    AbaqusPart part;
    part.name = params["NAME"];

    if (part.name.empty())
    {
        m_lastError = "Part name is required";
        return false;
    }

    feLogEx( nullptr, "  Parsing part: %s\n", part.name.c_str());

    // Parse part contents
    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);

        if (line.empty() || (line[0] == '*' && line[1] == '*'))
            continue;

        if (line[0] == '*')
        {
            std::string keyword = ToUpper(line.substr(1));

            if (keyword.find("NODE") == 0)
            {
                if (!ParseNode(file, part))
                    return false;
            }
            else if (keyword.find("ELEMENT") == 0)
            {
                if (!ParseElement(file, part))
                    return false;
            }
            else if (keyword.find("NSET") == 0)
            {
                if (!ParseNset(file, part))
                    return false;
            }
            else if (keyword.find("ELSET") == 0)
            {
                if (!ParseElset(file, part))
                    return false;
            }
            else if (keyword.find("SURFACE") == 0)
            {
                if (!ParseSurface(file, part))
                    return false;
            }
            else if (keyword.find("END PART") == 0)
            {
                break;
            }
        }
    }

    if (ValidatePart(part))
    {
        m_parts.push_back(part);
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseNode(std::ifstream& file, AbaqusPart& part)
{
    std::string line;
    int nodeCount = 0;

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);

        if (line.empty())
            continue;
        if (line[0] == '*')
        {
            // Put the line back for next keyword
            file.seekg(-(line.length() + 1), std::ios::cur);
            m_lineNumber--;
            break;
        }

        // Parse node data: ID, X, Y, Z
        std::vector<std::string> tokens = SplitString(line, ',');
        if (tokens.size() < 4)
        {
            m_lastError = "Invalid node data format";
            return false;
        }

        AbaqusNode node;
        node.id = std::stoi(TrimString(tokens[0]));
        node.x = std::stod(TrimString(tokens[1]));
        node.y = std::stod(TrimString(tokens[2]));
        node.z = std::stod(TrimString(tokens[3]));

        part.nodes.push_back(node);
        nodeCount++;
    }

    //feLog("    Read %d nodes\n", nodeCount);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseElement(std::ifstream& file, AbaqusPart& part)
{
    // Get element type from keyword line
    std::string keywordLine;
    file.seekg(-(1), std::ios::cur);
    std::getline(file, keywordLine);
    m_lineNumber++;

    std::map<std::string, std::string> params;
    ParseKeywordParams(keywordLine, params);

    std::string elemType = params["TYPE"];
    std::string elsetName = params["ELSET"];

    if (elemType.empty())
    {
        m_lastError = "Element type is required";
        return false;
    }

    int typeId = ConvertElementType(elemType);
    int nodeCount = GetElementNodeCount(typeId);
    int elemCount = 0;

    std::string line;
    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);

        if (line.empty())
            continue;
        if (line[0] == '*')
        {
            file.seekg(-(line.length() + 1), std::ios::cur);
            m_lineNumber--;
            break;
        }

        // Parse element connectivity
        std::vector<std::string> tokens = SplitString(line, ',');

        AbaqusElement elem;
        elem.id = std::stoi(TrimString(tokens[0]));
        elem.type = typeId;
        elem.elset = elsetName;

        // Read node connectivity
        for (size_t i = 1; i < tokens.size(); i++)
        {
            elem.nodes.push_back(std::stoi(TrimString(tokens[i])));
        }

        // Check if element spans multiple lines
        while (elem.nodes.size() < nodeCount && std::getline(file, line))
        {
            m_lineNumber++;
            line = TrimString(line);
            if (line.empty())
                continue;
            if (line[0] == '*')
            {
                file.seekg(-(line.length() + 1), std::ios::cur);
                m_lineNumber--;
                break;
            }

            tokens = SplitString(line, ',');
            for (const auto& token : tokens)
            {
                elem.nodes.push_back(std::stoi(TrimString(token)));
            }
        }

        part.elements.push_back(elem);
        elemCount++;
    }

    //feLog("    Read %d %s elements\n", elemCount, elemType.c_str());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseNset(std::ifstream& file, AbaqusPart& part)
{
    std::string keywordLine;
    file.seekg(-(1), std::ios::cur);
    std::getline(file, keywordLine);
    m_lineNumber++;

    std::map<std::string, std::string> params;
    ParseKeywordParams(keywordLine, params);

    std::string nsetName = params["NSET"];
    if (nsetName.empty())
    {
        m_lastError = "Node set name is required";
        return false;
    }

    std::vector<int> nodes;
    std::string line;

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);

        if (line.empty())
            continue;
        if (line[0] == '*')
        {
            file.seekg(-(line.length() + 1), std::ios::cur);
            m_lineNumber--;
            break;
        }

        std::vector<std::string> tokens = SplitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = TrimString(token);
            if (!trimmed.empty())
                nodes.push_back(std::stoi(trimmed));
        }
    }

    part.nodeSets[nsetName] = nodes;
    //feLog("    Created node set '%s' with %d nodes\n", nsetName.c_str(), (int)nodes.size());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseElset(std::ifstream& file, AbaqusPart& part)
{
    std::string keywordLine;
    file.seekg(-(1), std::ios::cur);
    std::getline(file, keywordLine);
    m_lineNumber++;

    std::map<std::string, std::string> params;
    ParseKeywordParams(keywordLine, params);

    std::string elsetName = params["ELSET"];
    if (elsetName.empty())
    {
        m_lastError = "Element set name is required";
        return false;
    }

    std::vector<int> elems;
    std::string line;

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);

        if (line.empty())
            continue;
        if (line[0] == '*')
        {
            file.seekg(-(line.length() + 1), std::ios::cur);
            m_lineNumber--;
            break;
        }

        std::vector<std::string> tokens = SplitString(line, ',');
        for (const auto& token : tokens)
        {
            std::string trimmed = TrimString(token);
            if (!trimmed.empty())
                elems.push_back(std::stoi(trimmed));
        }
    }

    part.elementSets[elsetName] = elems;
    //feLog("    Created element set '%s' with %d elements\n", elsetName.c_str(), (int)elems.size());
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseSurface(std::ifstream& file, AbaqusPart& part)
{
    // TODO: Implement surface parsing
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseMaterial(std::ifstream& file)
{
    // TODO: Implement material parsing
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseStep(std::ifstream& file)
{
    // TODO: Implement step parsing
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseAssembly(std::ifstream& file)
{
    // TODO: Implement assembly parsing for instances
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ProcessData(FEModel* fem)
{
    FEMesh& mesh = fem->GetMesh();

    // Count total nodes
    m_totalNodes = 0;
    for (const auto& part : m_parts)
    {
        m_totalNodes += part.nodes.size();
    }

    // Create nodes in mesh
    mesh.CreateNodes(m_totalNodes);

    // Process each part
    int nodeOffset = 0;
    for (const auto& part : m_parts)
    {
        //feLog("\nProcessing part: %s\n", part.name.c_str());

        // Create node mapping for this part
        std::map<int, int> nodeMap;
        if (!CreateNodes(part, &mesh, nodeMap))
            return false;

        // Create domain for this part
        RgDomain* domain = CreateDomain(part, &mesh);
        if (!domain)
            return false;

        mesh.AddDomain(domain);

        // Create elements
        if (!CreateElements(part, domain, nodeMap))
            return false;

        // Create node sets
        if (!CreateNodeSets(part, &mesh, nodeMap))
            return false;

        // Create element sets
        if (!CreateElementSets(part, &mesh, domain))
            return false;

        m_totalElements += part.elements.size();
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::CreateNodes(const AbaqusPart& part, FEMesh* mesh, std::map<int, int>& nodeMap)
{
    static int globalNodeIndex = 0;

    for (const auto& abqNode : part.nodes)
    {
        FENode& node = mesh->Node(globalNodeIndex);
        node.m_r0 = Vector3d(abqNode.x, abqNode.y, abqNode.z);
        node.m_rt = node.m_r0;
        node.SetID(abqNode.id);

        nodeMap[abqNode.id] = globalNodeIndex;
        m_globalNodeMap[abqNode.id] = globalNodeIndex;

        globalNodeIndex++;
    }

    return true;
}

//-----------------------------------------------------------------------------
RgDomain* AbaqusImport::CreateDomain(const AbaqusPart& part, FEMesh* mesh)
{
    // Determine domain type from first element
    if (part.elements.empty())
        return nullptr;

    // For now, create a solid domain
    // TODO: Support other domain types
    RgSolidDomain* domain = nullptr; //new RgSolidDomain(mesh->GetFEModel());
    domain->SetName(part.name);

    return domain;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::CreateElements(const AbaqusPart& part, RgDomain* domain, const std::map<int, int>& nodeMap)
{
    // TODO: Implement element creation
    // This requires knowledge of your element classes
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::CreateNodeSets(const AbaqusPart& part, FEMesh* mesh, const std::map<int, int>& nodeMap)
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

        // TODO: Set node indices to nodeSet
        mesh->AddNodeSet(nodeSet);
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::CreateElementSets(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain)
{
    for (const auto& elset : part.elementSets)
    {
        RgElementSet* elemSet = new RgElementSet(mesh->GetFEModel());
        elemSet->SetName(elset.first);

        // TODO: Add elements to set
        mesh->AddElementSet(elemSet);
    }

    return true;
}

//-----------------------------------------------------------------------------
int AbaqusImport::ConvertElementType(const std::string& abqType)
{
    std::string type = ToUpper(abqType);

    //// Solid elements
    //if (type == "C3D8" || type == "C3D8R")
    //    return ElementType::FE_HEX8G1;
    //if (type == "C3D20" || type == "C3D20R")
    //    return FE_HEX20;
    //if (type == "C3D27")
    //    return FE_HEX27;
    //if (type == "C3D4")
    //    return FE_TET4;
    //if (type == "C3D10")
    //    return FE_TET10;

    //// Shell elements
    //if (type == "S4" || type == "S4R")
    //    return FE_QUAD4;
    //if (type == "S8R")
    //    return FE_QUAD8;
    //if (type == "S3" || type == "S3R")
    //    return FE_TRI3;

    //// Beam/Truss
    //if (type == "T3D2")
    //    return FE_LINE2;

    return -1;  // Unknown type
}

//-----------------------------------------------------------------------------
int AbaqusImport::GetElementNodeCount(int elemType)
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
        case FE_QUAD4:
            return 4;
        case FE_QUAD8:
            return 8;
        case FE_TRI3:
            return 3;
        case FE_LINE2:
            return 2;*/
        default:
            return 0;
    }
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ValidatePart(const AbaqusPart& part)
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
std::string AbaqusImport::ReadKeywordLine(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);
    m_lineNumber++;
    return line;
}

//-----------------------------------------------------------------------------
void AbaqusImport::ParseKeywordParams(const std::string& line, std::map<std::string, std::string>& params)
{
    std::vector<std::string> tokens = SplitString(line, ',');

    for (size_t i = 1; i < tokens.size(); i++)
    {
        std::string token = TrimString(tokens[i]);
        size_t pos = token.find('=');

        if (pos != std::string::npos)
        {
            std::string key = ToUpper(TrimString(token.substr(0, pos)));
            std::string value = TrimString(token.substr(pos + 1));
            params[key] = value;
        }
    }
}

//-----------------------------------------------------------------------------
std::vector<std::string> AbaqusImport::SplitString(const std::string& str, char delimiter)
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
std::string AbaqusImport::TrimString(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

//-----------------------------------------------------------------------------
std::string AbaqusImport::ToUpper(const std::string& str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::SkipToNextKeyword(std::ifstream& file)
{
    std::string line;
    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = TrimString(line);
        if (!line.empty() && line[0] == '*')
        {
            file.seekg(-(line.length() + 1), std::ios::cur);
            m_lineNumber--;
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseHeading(std::ifstream& file)
{
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseInstance(std::ifstream& file)
{
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseBoundary(std::ifstream& file)
{
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::ParseLoad(std::ifstream& file)
{
    SkipToNextKeyword(file);
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusImport::CreateSurfaces(const AbaqusPart& part, FEMesh* mesh, RgDomain* domain)
{
    // TODO: Implement surface creation
    return true;
}