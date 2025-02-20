/*********************************************************************
 * \file   RgFemApp.h
 * \brief  
 * 
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "app/FEAppConfig.h"
#include "app/CmdOptions.h"

class FEModel;

class RgFemApp
{
public:
    static RgFemApp* Instance();

public:
    RgFemApp();

    bool Init(int nargs, char* argv[]);

    bool Configure(const char* szconfig);

    int Run();

    void Finish();

    bool ParseCmdLine(int nargs, char* argv[]);
    void ProcessCommands();

    // run an febio model
    int RunModel();

public:
    // get the current model
    FEBioModel* GetCurrentModel();

protected:
    // show FEBio prompt
    int prompt();

    // set the currently active model
    void SetCurrentModel(FEBioModel* fem);

    // apply configuration changes to model
    void ApplyConfig(FEBioModel& fem);

private:
    FEAppConfig m_config;  // configuration options
    CmdOptions m_cmd_opts;
    FEModel* mp_fem;
};
