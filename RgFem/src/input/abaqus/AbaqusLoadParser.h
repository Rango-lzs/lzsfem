/*********************************************************************
 * \file   AbaqusLoadParser.h
 * \brief  Parser for Abaqus load keywords (*Cload, *Dload)
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#pragma once
#include <string>
#include <vector>
#include <sstream>

class FEModel;
class RgLoad;
class FENodeSet;
class FEFacetSet;

//-----------------------------------------------------------------------------
//! Abaqus load parser
class AbaqusLoadParser
{
public:
    //! Load data structure
    struct LoadData
    {
        std::string name;           //!< Load name
        std::string type;           //!< Load type
        std::string loadType;       //!< Abaqus load type (P, GRAV, CENTRIF, etc)
        std::string targetName;     //!< Node set or element set name
        int dof;                    //!< DOF (for concentrated loads)
        double magnitude;           //!< Load magnitude
        double x, y, z;            //!< Components (for vector loads)
        bool isConcentrated;        //!< Is concentrated (nodal) load
        bool isDistributed;         //!< Is distributed (surface/body) load
        
        LoadData()
            : dof(0), magnitude(0.0), x(0.0), y(0.0), z(0.0),
              isConcentrated(false), isDistributed(false)
        {}
    };

public:
    AbaqusLoadParser(FEModel* fem);
    
    //! Parse *Cload (concentrated load) keyword
    //! \param lines Input file lines
    //! \param startIdx Starting line index (at *Cload keyword)
    //! \return Number of lines consumed
    int ParseCload(const std::vector<std::string>& lines, int startIdx);
    
    //! Parse *Dload (distributed load) keyword
    //! \param lines Input file lines
    //! \param startIdx Starting line index (at *Dload keyword)
    //! \return Number of lines consumed
    int ParseDload(const std::vector<std::string>& lines, int startIdx);
    
    //! Create loads from parsed data
    bool CreateLoads();
    
    //! Get parsed load data
    const std::vector<LoadData>& GetLoadData() const { return m_loadData; }

private:
    //! Parse *Cload keyword line
    bool ParseCloadKeyword(const std::string& line, std::string& loadName);
    
    //! Parse *Dload keyword line
    bool ParseDloadKeyword(const std::string& line, std::string& loadName);
    
    //! Parse *Cload data line
    bool ParseCloadData(const std::string& line, LoadData& data);
    
    //! Parse *Dload data line
    bool ParseDloadData(const std::string& line, LoadData& data);
    
    //! Convert Abaqus DOF to internal DOF
    int ConvertDOF(int abaqusDOF);
    
    //! Create load from parsed data
    RgLoad* CreateLoad(const LoadData& data);
    
    //! Trim whitespace
    std::string Trim(const std::string& str);
    
    //! Convert to uppercase
    std::string ToUpper(const std::string& str);
    
    //! Check if line is a keyword
    bool IsKeyword(const std::string& line);
    
    //! Parse vector from string
    bool ParseVector(const std::string& str, double& x, double& y, double& z);

private:
    FEModel* m_fem;
    std::vector<LoadData> m_loadData;
};

//-----------------------------------------------------------------------------
//! Abaqus load formats:
//!
//! *Cload (Concentrated Load - Nodal Forces):
//! ** Name: Load-1 Type: Concentrated force
//! *Cload
//! NodeSetName, DOF, Magnitude
//! 
//! Example:
//! *Cload
//! LoadPoint, 1, 100.0     # 100N in X direction
//! LoadPoint, 2, -50.0     # -50N in Y direction
//! 
//! DOF: 1,2,3 = X,Y,Z forces; 4,5,6 = X,Y,Z moments
//!
//! *Dload (Distributed Load - Pressure/Traction):
//! ** Name: Load-2 Type: Pressure
//! *Dload
//! SurfaceName, P, Magnitude
//! 
//! or with traction:
//! *Dload
//! SurfaceName, TRVEC, Magnitude, x, y, z
//! 
//! Load types:
//!   P      - Pressure (normal to surface)
//!   TRVEC  - Traction vector
//!   GRAV   - Gravity (x, y, z components)
//!   CENTRIF- Centrifugal (origin x,y,z, axis x,y,z, omega)
//!
//! Example pressure:
//! *Dload
//! TopFace, P, 1.0e6       # 1 MPa pressure
//!
//! Example gravity:
//! *Dload
//! AllElements, GRAV, 9.81, 0.0, 0.0, -1.0  # g in -Z direction
//!
//! Example centrifugal:
//! *Dload
//! AllElements, CENTRIF, 100.0, 0,0,0, 0,0,1  # 100 rad/s about Z axis
