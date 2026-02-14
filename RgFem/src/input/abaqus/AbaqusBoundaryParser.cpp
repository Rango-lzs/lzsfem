#include "AbaqusBoundaryParser.h"
#include "femcore/BoundaryCondition/RgBoundaryCondition.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENodeSet.h"
#include "logger/log.h"
#include <cctype>
#include <sstream>

//-----------------------------------------------------------------------------
AbaqusBoundaryParser::AbaqusBoundaryParser(FEModel* fem) 
    : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
int AbaqusBoundaryParser::ParseBoundary(const std::vector<std::string>& lines, int startIdx)
{
    if (startIdx >= lines.size()) return 0;
    
    std::string line = Trim(lines[startIdx]);
    
    // Should start with *Boundary
    if (ToUpper(line).find("*BOUNDARY") != 0)
    {
        RgLogError("Expected *Boundary keyword at line %d", startIdx + 1);
        return 0;
    }
    
    // Parse keyword line for name and type
    std::string bcName, bcType;
    ParseKeywordLine(line, bcName, bcType);
    
    // Look for comment lines with Name: and Type:
    // These might be before the *Boundary line
    if (startIdx > 0)
    {
        for (int i = startIdx - 1; i >= 0 && i >= startIdx - 5; --i)
        {
            std::string commentLine = Trim(lines[i]);
            if (commentLine.empty() || commentLine[0] != '*') continue;
            if (commentLine.size() < 2 || commentLine[1] != '*') break;
            
            // Parse comment
            std::string comment = Trim(commentLine.substr(2));
            
            // Look for "Name:"
            size_t namePos = comment.find("Name:");
            if (namePos != std::string::npos)
            {
                size_t typePos = comment.find("Type:", namePos);
                if (typePos != std::string::npos)
                {
                    bcName = Trim(comment.substr(namePos + 5, typePos - namePos - 5));
                    bcType = Trim(comment.substr(typePos + 5));
                }
                else
                {
                    bcName = Trim(comment.substr(namePos + 5));
                }
            }
        }
    }
    
    // Parse data lines
    int lineCount = 1;
    int currentLine = startIdx + 1;
    
    while (currentLine < lines.size())
    {
        line = Trim(lines[currentLine]);
        
        // Stop at empty line or next keyword
        if (line.empty() || IsKeyword(line))
            break;
        
        // Skip comments
        if (line[0] == '*' && line.size() > 1 && line[1] == '*')
        {
            currentLine++;
            lineCount++;
            continue;
        }
        
        // Parse boundary data
        BoundaryData data;
        data.name = bcName;
        data.type = bcType;
        
        if (ParseDataLine(line, data))
        {
            m_boundaryData.push_back(data);
        }
        else
        {
            RgLogWarning("Failed to parse boundary data at line %d: %s", 
                        currentLine + 1, line.c_str());
        }
        
        currentLine++;
        lineCount++;
    }
    
    return lineCount;
}

//-----------------------------------------------------------------------------
bool AbaqusBoundaryParser::ParseKeywordLine(const std::string& line, 
                                           std::string& bcName, 
                                           std::string& bcType)
{
    // *Boundary might have parameters like *Boundary, OP=NEW
    // We'll ignore parameters for now
    
    bcName = "";
    bcType = "";
    
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusBoundaryParser::ParseDataLine(const std::string& line, BoundaryData& data)
{
    // Format: NodeSetName, FirstDOF, LastDOF, Magnitude
    // Or:     NodeSetName, ENCASTRE
    
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    
    while (std::getline(ss, token, ','))
    {
        token = Trim(token);
        if (!token.empty())
            tokens.push_back(token);
    }
    
    if (tokens.empty())
        return false;
    
    // First token is node set name
    data.nodeSetName = tokens[0];
    
    if (tokens.size() == 1)
    {
        RgLogError("Boundary condition missing parameters for node set %s", 
                  data.nodeSetName.c_str());
        return false;
    }
    
    // Check for ENCASTRE
    std::string secondToken = ToUpper(tokens[1]);
    if (secondToken == "ENCASTRE")
    {
        data.isEncastre = true;
        data.firstDOF = 1;
        data.lastDOF = 6;
        data.magnitude = 0.0;
        return true;
    }
    
    // Parse DOF specification
    try
    {
        if (tokens.size() >= 2)
        {
            data.firstDOF = std::stoi(tokens[1]);
        }
        
        if (tokens.size() >= 3)
        {
            data.lastDOF = std::stoi(tokens[2]);
        }
        else
        {
            data.lastDOF = data.firstDOF;
        }
        
        if (tokens.size() >= 4)
        {
            data.magnitude = std::stod(tokens[3]);
        }
        else
        {
            data.magnitude = 0.0;  // Fixed BC
        }
    }
    catch (const std::exception& e)
    {
        RgLogError("Error parsing boundary data: %s", e.what());
        return false;
    }
    
    // Validate DOF range
    if (data.firstDOF < 1 || data.firstDOF > 6)
    {
        RgLogError("Invalid first DOF: %d (must be 1-6)", data.firstDOF);
        return false;
    }
    
    if (data.lastDOF < data.firstDOF || data.lastDOF > 6)
    {
        RgLogError("Invalid last DOF: %d", data.lastDOF);
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusBoundaryParser::CreateBoundaryConditions()
{
    for (size_t i = 0; i < m_boundaryData.size(); ++i)
    {
        const BoundaryData& data = m_boundaryData[i];
        
        // Find node set
        FENodeSet* nodeSet = m_fem->GetMesh().FindNodeSet(data.nodeSetName);
        if (!nodeSet)
        {
            RgLogError("Node set '%s' not found for boundary condition", 
                      data.nodeSetName.c_str());
            return false;
        }
        
        // Create appropriate BC
        if (data.isEncastre)
        {
            // ENCASTRE: fix all DOFs
            RgFixedBC* bc = RgBCFactory::CreateEncastre(m_fem, nodeSet);
            bc->SetName(data.name.empty() ? data.nodeSetName + "_encastre" : data.name);
            m_fem->AddBoundaryCondition(bc);
            
            RgLog("Created ENCASTRE BC on node set '%s'\n", data.nodeSetName.c_str());
        }
        else
        {
            // Check if it's a fixed BC (magnitude = 0) or prescribed BC
            bool isFixed = (fabs(data.magnitude) < 1e-20);
            
            // Create BC for each DOF in range
            for (int dof = data.firstDOF; dof <= data.lastDOF; ++dof)
            {
                int internalDOF = ConvertDOF(dof);
                
                if (isFixed)
                {
                    // Fixed BC
                    int dofMask = (1 << internalDOF);
                    RgFixedBC* bc = RgBCFactory::CreateFixed(m_fem, nodeSet, dofMask);
                    
                    std::string name = data.name.empty() ? 
                        data.nodeSetName + "_fixed_" + std::to_string(dof) : 
                        data.name + "_" + std::to_string(dof);
                    bc->SetName(name);
                    
                    m_fem->AddBoundaryCondition(bc);
                    
                    RgLog("Created fixed BC on node set '%s', DOF %d\n", 
                         data.nodeSetName.c_str(), dof);
                }
                else
                {
                    // Prescribed BC
                    if (internalDOF < 3)
                    {
                        // Displacement
                        RgPrescribedDisplacement* bc = RgBCFactory::CreatePrescribed(
                            m_fem, nodeSet, internalDOF, data.magnitude);
                        
                        std::string name = data.name.empty() ? 
                            data.nodeSetName + "_disp_" + std::to_string(dof) : 
                            data.name + "_" + std::to_string(dof);
                        bc->SetName(name);
                        
                        m_fem->AddBoundaryCondition(bc);
                        
                        RgLog("Created prescribed displacement BC on node set '%s', "
                             "DOF %d, magnitude %g\n", 
                             data.nodeSetName.c_str(), dof, data.magnitude);
                    }
                    else
                    {
                        // Rotation
                        RgPrescribedRotation* bc = new RgPrescribedRotation(m_fem);
                        bc->SetNodeSet(nodeSet);
                        bc->SetDOF(internalDOF);
                        bc->SetScale(data.magnitude);
                        
                        std::string name = data.name.empty() ? 
                            data.nodeSetName + "_rot_" + std::to_string(dof) : 
                            data.name + "_" + std::to_string(dof);
                        bc->SetName(name);
                        
                        m_fem->AddBoundaryCondition(bc);
                        
                        RgLog("Created prescribed rotation BC on node set '%s', "
                             "DOF %d, magnitude %g\n", 
                             data.nodeSetName.c_str(), dof, data.magnitude);
                    }
                }
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
int AbaqusBoundaryParser::ConvertDOF(int abaqusDOF)
{
    // Abaqus: 1,2,3 = X,Y,Z displacement; 4,5,6 = X,Y,Z rotation
    // Internal: 0,1,2 = X,Y,Z displacement; 3,4,5 = X,Y,Z rotation
    
    if (abaqusDOF < 1 || abaqusDOF > 6)
    {
        RgLogWarning("Invalid Abaqus DOF: %d, using 0", abaqusDOF);
        return 0;
    }
    
    return abaqusDOF - 1;
}

//-----------------------------------------------------------------------------
RgBoundaryCondition* AbaqusBoundaryParser::CreateBC(const BoundaryData& data)
{
    // This is a helper that could be extended for different BC types
    return nullptr;
}

//-----------------------------------------------------------------------------
std::string AbaqusBoundaryParser::Trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

//-----------------------------------------------------------------------------
std::string AbaqusBoundaryParser::ToUpper(const std::string& str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

//-----------------------------------------------------------------------------
bool AbaqusBoundaryParser::IsKeyword(const std::string& line)
{
    if (line.empty()) return false;
    
    std::string trimmed = Trim(line);
    if (trimmed.empty()) return false;
    
    // Keyword lines start with '*' but not '**' (comment)
    if (trimmed[0] == '*')
    {
        if (trimmed.size() > 1 && trimmed[1] == '*')
            return false;  // Comment
        return true;
    }
    
    return false;
}
