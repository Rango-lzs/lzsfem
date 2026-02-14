#include "AbaqusLoadParser.h"
#include "femcore/Load/RgLoad.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include "femcore/FENodeSet.h"
#include "femcore/FEFacetSet.h"
#include "femcore/FESurface.h"
#include "logger/log.h"
#include <cctype>
#include <sstream>
#include <algorithm>

//-----------------------------------------------------------------------------
AbaqusLoadParser::AbaqusLoadParser(FEModel* fem)
    : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
int AbaqusLoadParser::ParseCload(const std::vector<std::string>& lines, int startIdx)
{
    if (startIdx >= lines.size()) return 0;
    
    std::string line = Trim(lines[startIdx]);
    
    // Should start with *Cload
    if (ToUpper(line).find("*CLOAD") != 0)
    {
        RgLogError("Expected *Cload keyword at line %d", startIdx + 1);
        return 0;
    }
    
    // Parse keyword line
    std::string loadName;
    ParseCloadKeyword(line, loadName);
    
    // Look for comment lines with Name: and Type:
    if (startIdx > 0)
    {
        for (int i = startIdx - 1; i >= 0 && i >= startIdx - 5; --i)
        {
            std::string commentLine = Trim(lines[i]);
            if (commentLine.empty() || commentLine[0] != '*') continue;
            if (commentLine.size() < 2 || commentLine[1] != '*') break;
            
            std::string comment = Trim(commentLine.substr(2));
            
            size_t namePos = comment.find("Name:");
            if (namePos != std::string::npos)
            {
                size_t typePos = comment.find("Type:", namePos);
                if (typePos != std::string::npos)
                {
                    loadName = Trim(comment.substr(namePos + 5, typePos - namePos - 5));
                }
                else
                {
                    loadName = Trim(comment.substr(namePos + 5));
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
        
        if (line.empty() || IsKeyword(line))
            break;
        
        if (line[0] == '*' && line.size() > 1 && line[1] == '*')
        {
            currentLine++;
            lineCount++;
            continue;
        }
        
        LoadData data;
        data.name = loadName;
        data.type = "Concentrated";
        data.isConcentrated = true;
        
        if (ParseCloadData(line, data))
        {
            m_loadData.push_back(data);
        }
        else
        {
            RgLogWarning("Failed to parse Cload data at line %d: %s",
                        currentLine + 1, line.c_str());
        }
        
        currentLine++;
        lineCount++;
    }
    
    return lineCount;
}

//-----------------------------------------------------------------------------
int AbaqusLoadParser::ParseDload(const std::vector<std::string>& lines, int startIdx)
{
    if (startIdx >= lines.size()) return 0;
    
    std::string line = Trim(lines[startIdx]);
    
    if (ToUpper(line).find("*DLOAD") != 0)
    {
        RgLogError("Expected *Dload keyword at line %d", startIdx + 1);
        return 0;
    }
    
    std::string loadName;
    ParseDloadKeyword(line, loadName);
    
    // Look for comment lines
    if (startIdx > 0)
    {
        for (int i = startIdx - 1; i >= 0 && i >= startIdx - 5; --i)
        {
            std::string commentLine = Trim(lines[i]);
            if (commentLine.empty() || commentLine[0] != '*') continue;
            if (commentLine.size() < 2 || commentLine[1] != '*') break;
            
            std::string comment = Trim(commentLine.substr(2));
            
            size_t namePos = comment.find("Name:");
            if (namePos != std::string::npos)
            {
                size_t typePos = comment.find("Type:", namePos);
                if (typePos != std::string::npos)
                {
                    loadName = Trim(comment.substr(namePos + 5, typePos - namePos - 5));
                }
                else
                {
                    loadName = Trim(comment.substr(namePos + 5));
                }
            }
        }
    }
    
    int lineCount = 1;
    int currentLine = startIdx + 1;
    
    while (currentLine < lines.size())
    {
        line = Trim(lines[currentLine]);
        
        if (line.empty() || IsKeyword(line))
            break;
        
        if (line[0] == '*' && line.size() > 1 && line[1] == '*')
        {
            currentLine++;
            lineCount++;
            continue;
        }
        
        LoadData data;
        data.name = loadName;
        data.type = "Distributed";
        data.isDistributed = true;
        
        if (ParseDloadData(line, data))
        {
            m_loadData.push_back(data);
        }
        else
        {
            RgLogWarning("Failed to parse Dload data at line %d: %s",
                        currentLine + 1, line.c_str());
        }
        
        currentLine++;
        lineCount++;
    }
    
    return lineCount;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::ParseCloadKeyword(const std::string& line, std::string& loadName)
{
    loadName = "";
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::ParseDloadKeyword(const std::string& line, std::string& loadName)
{
    loadName = "";
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::ParseCloadData(const std::string& line, LoadData& data)
{
    // Format: NodeSetName, DOF, Magnitude
    
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    
    while (std::getline(ss, token, ','))
    {
        token = Trim(token);
        if (!token.empty())
            tokens.push_back(token);
    }
    
    if (tokens.size() < 3)
    {
        RgLogError("Cload data requires at least 3 fields: NodeSet, DOF, Magnitude");
        return false;
    }
    
    try
    {
        data.targetName = tokens[0];
        data.dof = std::stoi(tokens[1]);
        data.magnitude = std::stod(tokens[2]);
    }
    catch (const std::exception& e)
    {
        RgLogError("Error parsing Cload data: %s", e.what());
        return false;
    }
    
    // Validate DOF
    if (data.dof < 1 || data.dof > 6)
    {
        RgLogError("Invalid DOF: %d (must be 1-6)", data.dof);
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::ParseDloadData(const std::string& line, LoadData& data)
{
    // Format: SurfaceName, LoadType, Magnitude [, components...]
    
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    
    while (std::getline(ss, token, ','))
    {
        token = Trim(token);
        if (!token.empty())
            tokens.push_back(token);
    }
    
    if (tokens.size() < 3)
    {
        RgLogError("Dload data requires at least 3 fields");
        return false;
    }
    
    data.targetName = tokens[0];
    data.loadType = ToUpper(tokens[1]);
    
    try
    {
        data.magnitude = std::stod(tokens[2]);
    }
    catch (const std::exception& e)
    {
        RgLogError("Error parsing magnitude: %s", e.what());
        return false;
    }
    
    // Parse additional components based on load type
    if (data.loadType == "P")
    {
        // Pressure - no additional components
    }
    else if (data.loadType == "TRVEC")
    {
        // Traction vector - needs x, y, z components
        if (tokens.size() >= 6)
        {
            data.x = std::stod(tokens[3]);
            data.y = std::stod(tokens[4]);
            data.z = std::stod(tokens[5]);
        }
        else
        {
            RgLogError("TRVEC requires x, y, z components");
            return false;
        }
    }
    else if (data.loadType == "GRAV")
    {
        // Gravity - magnitude is acceleration, x,y,z is direction
        if (tokens.size() >= 6)
        {
            data.x = std::stod(tokens[3]);
            data.y = std::stod(tokens[4]);
            data.z = std::stod(tokens[5]);
        }
        else
        {
            RgLogError("GRAV requires x, y, z direction components");
            return false;
        }
    }
    else if (data.loadType == "CENTRIF")
    {
        // Centrifugal - magnitude is omega^2, then origin and axis
        // Format: name, CENTRIF, omega^2, x0, y0, z0, x1, y1, z1
        // where (x0,y0,z0) is origin and (x1,y1,z1) is axis direction
        if (tokens.size() >= 9)
        {
            // Origin stored in x,y,z
            data.x = std::stod(tokens[3]);
            data.y = std::stod(tokens[4]);
            data.z = std::stod(tokens[5]);
            
            // Axis direction stored separately (would need extension)
            // For now, store in load-specific manner
        }
        else
        {
            RgLogError("CENTRIF requires origin (x,y,z) and axis (x,y,z)");
            return false;
        }
    }
    else
    {
        RgLogWarning("Unknown Dload type: %s", data.loadType.c_str());
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::CreateLoads()
{
    for (size_t i = 0; i < m_loadData.size(); ++i)
    {
        const LoadData& data = m_loadData[i];
        
        RgLoad* load = CreateLoad(data);
        if (!load)
        {
            RgLogError("Failed to create load from data: %s", data.name.c_str());
            return false;
        }
        
        // Set name
        std::string name = data.name.empty() ? 
            data.targetName + "_load" : data.name;
        load->SetName(name);
        
        // Add to model
        m_fem->AddModelLoad(load);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
RgLoad* AbaqusLoadParser::CreateLoad(const LoadData& data)
{
    FEMesh& mesh = m_fem->GetMesh();
    
    if (data.isConcentrated)
    {
        // Concentrated (nodal) load
        FENodeSet* nodeSet = mesh.FindNodeSet(data.targetName);
        if (!nodeSet)
        {
            RgLogError("Node set '%s' not found for Cload", data.targetName.c_str());
            return nullptr;
        }
        
        int internalDOF = ConvertDOF(data.dof);
        
        if (internalDOF < 3)
        {
            // Force
            RgNodalLoad* load = RgLoadFactory::CreateNodalLoad(
                m_fem, nodeSet, internalDOF, data.magnitude);
            
            RgLog("Created nodal force load on '%s', DOF %d, magnitude %g\n",
                 data.targetName.c_str(), data.dof, data.magnitude);
            
            return load;
        }
        else
        {
            // Moment
            RgMomentLoad* load = new RgMomentLoad(m_fem);
            load->SetNodeSet(nodeSet);
            
            Vector3d moment(0, 0, 0);
            moment(internalDOF - 3) = data.magnitude;
            load->SetMoment(moment);
            load->SetMagnitude(1.0);
            
            RgLog("Created moment load on '%s', DOF %d, magnitude %g\n",
                 data.targetName.c_str(), data.dof, data.magnitude);
            
            return load;
        }
    }
    else if (data.isDistributed)
    {
        // Distributed load
        if (data.loadType == "P")
        {
            // Pressure load
            FEFacetSet* facetSet = mesh.FindFacetSet(data.targetName);
            if (!facetSet)
            {
                RgLogError("Facet set '%s' not found for pressure load", 
                          data.targetName.c_str());
                return nullptr;
            }
            
            RgSurfaceLoad* load = RgLoadFactory::CreatePressure(
                m_fem, facetSet, data.magnitude);
            
            RgLog("Created pressure load on '%s', magnitude %g\n",
                 data.targetName.c_str(), data.magnitude);
            
            return load;
        }
        else if (data.loadType == "TRVEC")
        {
            // Traction vector
            FEFacetSet* facetSet = mesh.FindFacetSet(data.targetName);
            if (!facetSet)
            {
                RgLogError("Facet set '%s' not found for traction load",
                          data.targetName.c_str());
                return nullptr;
            }
            
            FESurface* surface = mesh.CreateSurface(*facetSet);
            Vector3d traction(data.x, data.y, data.z);
            
            RgSurfaceLoad* load = RgLoadFactory::CreateTraction(
                m_fem, surface, traction);
            load->SetMagnitude(data.magnitude);
            
            RgLog("Created traction load on '%s', traction (%g, %g, %g), magnitude %g\n",
                 data.targetName.c_str(), data.x, data.y, data.z, data.magnitude);
            
            return load;
        }
        else if (data.loadType == "GRAV")
        {
            // Gravity load
            Vector3d g(data.x, data.y, data.z);
            g.unit();  // Normalize direction
            g *= data.magnitude;  // Scale by acceleration
            
            RgBodyLoad* load = RgLoadFactory::CreateGravity(m_fem, g);
            
            RgLog("Created gravity load, g = (%g, %g, %g) m/s^2\n",
                 g.x, g.y, g.z);
            
            return load;
        }
        else if (data.loadType == "CENTRIF")
        {
            // Centrifugal load
            Vector3d origin(data.x, data.y, data.z);
            // Note: Axis direction would need to be stored in extended LoadData
            // For now, assume Z-axis
            Vector3d axis(0, 0, 1);
            
            // Abaqus stores omega^2, we need omega
            double omega = sqrt(data.magnitude);
            
            RgBodyLoad* load = RgLoadFactory::CreateCentrifugal(
                m_fem, axis, origin, omega);
            
            RgLog("Created centrifugal load, omega = %g rad/s, origin = (%g, %g, %g)\n",
                 omega, origin.x, origin.y, origin.z);
            
            return load;
        }
        else
        {
            RgLogError("Unsupported Dload type: %s", data.loadType.c_str());
            return nullptr;
        }
    }
    
    return nullptr;
}

//-----------------------------------------------------------------------------
int AbaqusLoadParser::ConvertDOF(int abaqusDOF)
{
    // Abaqus: 1,2,3 = X,Y,Z forces; 4,5,6 = X,Y,Z moments
    // Internal: 0,1,2 = X,Y,Z forces; 3,4,5 = X,Y,Z moments
    
    if (abaqusDOF < 1 || abaqusDOF > 6)
    {
        RgLogWarning("Invalid Abaqus DOF: %d, using 0", abaqusDOF);
        return 0;
    }
    
    return abaqusDOF - 1;
}

//-----------------------------------------------------------------------------
std::string AbaqusLoadParser::Trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

//-----------------------------------------------------------------------------
std::string AbaqusLoadParser::ToUpper(const std::string& str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::IsKeyword(const std::string& line)
{
    if (line.empty()) return false;
    
    std::string trimmed = Trim(line);
    if (trimmed.empty()) return false;
    
    if (trimmed[0] == '*')
    {
        if (trimmed.size() > 1 && trimmed[1] == '*')
            return false;  // Comment
        return true;
    }
    
    return false;
}

//-----------------------------------------------------------------------------
bool AbaqusLoadParser::ParseVector(const std::string& str, double& x, double& y, double& z)
{
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, ','))
    {
        token = Trim(token);
        if (!token.empty())
            tokens.push_back(token);
    }
    
    if (tokens.size() < 3)
        return false;
    
    try
    {
        x = std::stod(tokens[0]);
        y = std::stod(tokens[1]);
        z = std::stod(tokens[2]);
    }
    catch (const std::exception&)
    {
        return false;
    }
    
    return true;
}
