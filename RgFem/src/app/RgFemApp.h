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
    int Run();
    void Finish();

    bool ParseCmdLine(int nargs, char* argv[]);
    int RunModel();

    FEModel* GetCurrentModel();
    void SetCurrentModel(FEModel* fem);

private:
    FEAppConfig m_config;  // configuration options
    CmdOptions m_cmd_opts;
    FEModel* mp_model;
};


#define GET_FEMODEL RgFemApp::Instance()->GetCurrentModel()