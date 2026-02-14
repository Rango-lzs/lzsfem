/*********************************************************************
 * \file   AbaqusBoundaryParser.h
 * \brief  Parser for Abaqus *Boundary keyword
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

class FEModel;
class RgBoundaryCondition;
class FENodeSet;

//-----------------------------------------------------------------------------
//! Abaqus boundary condition parser
class AbaqusBoundaryParser
{
public:
    //! Parse result structure
    struct BoundaryData
    {
        std::string name;           //!< BC name
        std::string type;           //!< BC type (ENCASTRE, displacement, etc)
        std::string nodeSetName;    //!< Node set name
        int firstDOF;               //!< First DOF (1-based, Abaqus convention)
        int lastDOF;                //!< Last DOF (1-based)
        double magnitude;           //!< Prescribed value
        bool isEncastre;            //!< Is ENCASTRE type
        bool isSymmetry;            //!< Is symmetry/antisymmetry
        
        BoundaryData() 
            : firstDOF(0), lastDOF(0), magnitude(0.0), 
              isEncastre(false), isSymmetry(false) 
        {}
    };

public:
    AbaqusBoundaryParser(FEModel* fem);
    
    //! Parse boundary conditions from input lines
    //! \param lines Input file lines
    //! \param startIdx Starting line index (at *Boundary keyword)
    //! \return Number of lines consumed
    int ParseBoundary(const std::vector<std::string>& lines, int startIdx);
    
    //! Create FE boundary conditions from parsed data
    bool CreateBoundaryConditions();
    
    //! Get parsed boundary data
    const std::vector<BoundaryData>& GetBoundaryData() const { return m_boundaryData; }

private:
    //! Parse *Boundary keyword line
    bool ParseKeywordLine(const std::string& line, std::string& bcName, 
                         std::string& bcType);
    
    //! Parse boundary data line
    bool ParseDataLine(const std::string& line, BoundaryData& data);
    
    //! Convert Abaqus DOF to internal DOF (0-based)
    int ConvertDOF(int abaqusDOF);
    
    //! Create BC from parsed data
    RgBoundaryCondition* CreateBC(const BoundaryData& data);
    
    //! Trim whitespace from string
    std::string Trim(const std::string& str);
    
    //! Convert string to uppercase
    std::string ToUpper(const std::string& str);
    
    //! Check if line is a keyword
    bool IsKeyword(const std::string& line);

private:
    FEModel* m_fem;
    std::vector<BoundaryData> m_boundaryData;
};

//-----------------------------------------------------------------------------
//! Abaqus *Boundary format:
//! 
//! ** BOUNDARY CONDITIONS **
//! ** Name: BC-1 Type: Symmetry/Antisymmetry/Encastre
//! *Boundary
//! Set-1, ENCASTRE
//! 
//! ** Name: BC-2 Type: Displacement/Rotation
//! *Boundary
//! Set-5, 1, 1, 0.1
//! Set-5, 2, 2, 0.1
//! Set-5, 3, 3, 0.1
//! 
//! Format: NodeSetName, FirstDOF, LastDOF, Magnitude
//! DOF numbering (Abaqus convention):
//!   1, 2, 3 = X, Y, Z displacement
//!   4, 5, 6 = X, Y, Z rotation
//! 
//! ENCASTRE: All DOFs fixed (displacement and rotation)
