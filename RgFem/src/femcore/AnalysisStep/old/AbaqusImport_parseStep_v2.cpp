/*********************************************************************
 * \file   AbaqusImport_parseStep_v2.cpp
 * \brief  Updated parseStep implementation with RgAnalysisStep support
 *
 * \author 
 * \date   February 2025
 *********************************************************************/

#include "AbaqusImport.h"
#include "femcore/FEModel.h"
#include "RgAnalysisStep.h"
#include "femcore/FEBoundaryCondition.h"
#include "femcore/FEModelLoad.h"
#include "logger/log.h"
#include <algorithm>
#include <sstream>

//-----------------------------------------------------------------------------
bool AbaqusImport::parseStep(std::ifstream& file, const std::string& keywordLine)
{
    RgLog("=================================================\n");
    RgLog("Parsing Step Section\n");
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

    // Create the analysis step object
    RgAnalysisStep* analysisStep = nullptr;

    // Temporary storage for BCs and loads defined in this step
    std::vector<FEBoundaryCondition*> stepBCs;
    std::vector<FEModelLoad*> stepLoads;

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
                RgLog("  Procedure: STATIC\n");

                // Create static analysis step
                if (!analysisStep)
                {
                    analysisStep = new RgStaticAnalysisStep(m_model);
                    analysisStep->SetName(stepInfo.name);
                    analysisStep->SetStepNumber(m_currentStep);
                }

                RgStaticAnalysisStep* staticStep = dynamic_cast<RgStaticAnalysisStep*>(analysisStep);
                if (staticStep)
                {
                    // Parse STATIC parameters
                    std::map<std::string, std::string> procParams;
                    parseKeywordParams(line, procParams);

                    // Check for NLGEOM (nonlinear geometry)
                    if (procParams.find("NLGEOM") != procParams.end())
                    {
                        std::string nlgeom = toUpper(procParams["NLGEOM"]);
                        staticStep->SetNonlinear(nlgeom == "YES" || nlgeom == "ON");
                    }
                }

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
                            RgLog("    Initial time increment: %g\n", stepInfo.initialTimeIncrement);
                        }
                        if (tokens.size() >= 2)
                        {
                            stepInfo.timePeriod = std::stod(trimString(tokens[1]));
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
                RgLog("  Procedure: DYNAMIC\n");

                // Create dynamic analysis step
                if (!analysisStep)
                {
                    analysisStep = new RgDynamicAnalysisStep(m_model);
                    analysisStep->SetName(stepInfo.name);
                    analysisStep->SetStepNumber(m_currentStep);
                }

                RgDynamicAnalysisStep* dynamicStep = dynamic_cast<RgDynamicAnalysisStep*>(analysisStep);
                if (dynamicStep)
                {
                    // Parse DYNAMIC parameters
                    std::map<std::string, std::string> procParams;
                    parseKeywordParams(line, procParams);

                    // Check for EXPLICIT
                    if (procParams.find("EXPLICIT") != procParams.end())
                    {
                        dynamicStep->SetExplicit(true);
                        RgLog("    Using explicit integration\n");
                    }
                }

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
                            RgLog("    Initial time increment: %g\n", stepInfo.initialTimeIncrement);
                        }
                        if (tokens.size() >= 2)
                        {
                            stepInfo.timePeriod = std::stod(trimString(tokens[1]));
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
                
                // Store current BC count
                int bcCountBefore = m_model->BoundaryConditions();
                
                // Parse boundary conditions
                if (!parseBoundary(file))
                {
                    RgLogWarning("    Failed to parse boundary conditions\n");
                }
                
                // Collect BCs added during parsing
                int bcCountAfter = m_model->BoundaryConditions();
                for (int i = bcCountBefore; i < bcCountAfter; ++i)
                {
                    FEBoundaryCondition* bc = m_model->BoundaryCondition(i);
                    stepBCs.push_back(bc);
                    RgLog("    Added BC '%s' to step\n", bc->GetName().c_str());
                }
            }
            else if (keyword.find("CLOAD") == 0)
            {
                RgLog("  Parsing CLOAD in step\n");
                
                // Store current load count
                int loadCountBefore = m_model->ModelLoads();
                
                // Parse concentrated loads
                if (!parseCload(file))
                {
                    RgLogWarning("    Failed to parse concentrated loads\n");
                }
                
                // Collect loads added during parsing
                int loadCountAfter = m_model->ModelLoads();
                for (int i = loadCountBefore; i < loadCountAfter; ++i)
                {
                    FEModelLoad* load = m_model->ModelLoad(i);
                    stepLoads.push_back(load);
                    RgLog("    Added load '%s' to step\n", load->GetName().c_str());
                }
            }
            else if (keyword.find("DLOAD") == 0)
            {
                RgLog("  Parsing DLOAD in step\n");
                
                // Store current load count
                int loadCountBefore = m_model->ModelLoads();
                
                // Parse distributed loads
                if (!parseDload(file))
                {
                    RgLogWarning("    Failed to parse distributed loads\n");
                }
                
                // Collect loads added during parsing
                int loadCountAfter = m_model->ModelLoads();
                for (int i = loadCountBefore; i < loadCountAfter; ++i)
                {
                    FEModelLoad* load = m_model->ModelLoad(i);
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

    // If no procedure was specified, create default static step
    if (!analysisStep)
    {
        analysisStep = new RgStaticAnalysisStep(m_model);
        analysisStep->SetName(stepInfo.name);
        analysisStep->SetStepNumber(m_currentStep);
        RgLog("  No procedure specified, using default STATIC\n");
    }

    // Set time parameters
    analysisStep->SetTimePeriod(stepInfo.timePeriod);
    analysisStep->SetInitialTimeIncrement(stepInfo.initialTimeIncrement);
    analysisStep->SetMinTimeIncrement(stepInfo.minTimeIncrement);
    analysisStep->SetMaxTimeIncrement(stepInfo.maxTimeIncrement);

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
        RgAnalysisStep* prevStep = dynamic_cast<RgAnalysisStep*>(m_model->GetStep(m_currentStep - 1));
        if (prevStep)
        {
            analysisStep->SetPreviousStep(prevStep);
            analysisStep->SetActivationMode(StepActivationMode::INHERITED);
            RgLog("  Inheriting BCs and loads from previous step\n");
        }
    }
    else
    {
        // First step - no inheritance
        analysisStep->SetActivationMode(StepActivationMode::NEW);
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
    RgLog("=================================================\n\n");

    return true;
}

//-----------------------------------------------------------------------------
// Helper function to skip to next keyword
bool AbaqusImport::skipToNextKeyword(std::ifstream& file)
{
    std::string line;
    std::streampos lastPos = file.tellg();

    while (std::getline(file, line))
    {
        m_lineNumber++;
        line = trimString(line);
        
        if (!line.empty() && line[0] == '*' && (line.size() < 2 || line[1] != '*'))
        {
            // Found next keyword, seek back
            file.seekg(lastPos);
            m_lineNumber--;
            return true;
        }
        
        lastPos = file.tellg();
    }
    
    return false;
}
